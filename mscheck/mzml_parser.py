"""
Minimal mzML parser for MSCheck.

Reads three data channels from mzML files using only stdlib + numpy:
  1. MS1 mass spectra   — m/z and intensity arrays, RT, polarity
  2. 2D DAD spectra     — full UV spectrum (wavelength vs absorbance) at every scan
  3. UV chromatograms   — single-wavelength intensity traces stored in <chromatogramList>

Replaces pyopenms, which could not read the 2D DAD electromagnetic radiation spectra
(MS:1000804) because it only decodes m/z arrays (MS:1000514), not wavelength arrays
(MS:1000617).
"""

from __future__ import annotations

import base64
import re
import zlib
import xml.etree.ElementTree as ET
from typing import Optional

import numpy as np

from .spectra_types import (
    ExperimentData,
    MS1Spectrum,
    UVChromatogram,
    UVSpectraMatrix,
)

# Backward-compatible alias
MzMLData = ExperimentData

# mzML namespace — may be absent in some files, handled via _find() / _iter()
_NS = "http://psi.hupo.org/ms/mzml"

# CV accessions referenced in parsing
_CV_MS1          = "MS:1000579"  # MS1 spectrum
_CV_EM_RAD       = "MS:1000804"  # electromagnetic radiation spectrum (DAD 2D)
_CV_POSITIVE     = "MS:1000130"  # positive scan
_CV_NEGATIVE     = "MS:1000129"  # negative scan
_CV_SCAN_TIME    = "MS:1000016"  # scan start time
_CV_MZ_ARRAY     = "MS:1000514"  # m/z array
_CV_INTENSITY    = "MS:1000515"  # intensity array
_CV_WAVELENGTH   = "MS:1000617"  # wavelength array
_CV_FLOAT64      = "MS:1000523"  # 64-bit float
_CV_FLOAT32      = "MS:1000521"  # 32-bit float
_CV_ZLIB         = "MS:1000574"  # zlib compression
_CV_NO_COMPRESS  = "MS:1000576"  # no compression


# ── Internal helpers ──────────────────────────────────────────────────────────

def _tag(local: str) -> str:
    """Return both namespaced and plain tag strings for find/iter calls."""
    return local  # handled via _find / _iter wrappers below


def _find(elem: ET.Element, local: str) -> Optional[ET.Element]:
    """Find a direct child by local tag name, with or without namespace."""
    # Explicit `is not None` — Element with no children is falsy on Py<=3.11.
    found = elem.find(f"{{{_NS}}}{local}")
    if found is not None:
        return found
    return elem.find(local)


def _iter(elem: ET.Element, local: str):
    """Iterate over direct children matching local tag name."""
    yield from elem.iterfind(f"{{{_NS}}}{local}")
    yield from elem.iterfind(local)


def _has_cv(elem: ET.Element, accession: str) -> bool:
    """Return True if *elem* has a direct cvParam child with *accession*."""
    for child in elem:
        if child.get("accession") == accession:
            return True
    return False


def _cv_value(elem: ET.Element, accession: str) -> Optional[str]:
    """Return the value attribute of a direct cvParam child, or None."""
    for child in elem:
        if child.get("accession") == accession:
            return child.get("value")
    return None


def _decode_binary_array(bda_elem: ET.Element) -> np.ndarray:
    """
    Decode one <binaryDataArray> element.

    Handles:
      dtype       : 64-bit float (default) or 32-bit float
      compression : zlib or none
    """
    is_32bit = _has_cv(bda_elem, _CV_FLOAT32)
    is_zlib  = _has_cv(bda_elem, _CV_ZLIB)

    binary_elem = _find(bda_elem, "binary")
    text = (binary_elem.text or "").strip() if binary_elem is not None else ""
    if not text:
        return np.array([], dtype=np.float64)

    raw = base64.b64decode(text)
    if is_zlib:
        raw = zlib.decompress(raw)

    dtype = np.float32 if is_32bit else np.float64
    return np.frombuffer(raw, dtype=dtype).copy()


def _extract_rt(spectrum_elem: ET.Element) -> float:
    """Extract scan start time from a <spectrum> element, returned in minutes."""
    scan_list = _find(spectrum_elem, "scanList")
    if scan_list is None:
        return 0.0
    for scan in _iter(scan_list, "scan"):
        val = _cv_value(scan, _CV_SCAN_TIME)
        if val is not None:
            # Check unit: UO:0000031 = minute, UO:0000010 = second
            rt = float(val)
            unit = None
            for child in scan:
                if child.get("accession") == _CV_SCAN_TIME:
                    unit = child.get("unitAccession", "")
                    break
            if unit == "UO:0000010":  # seconds
                rt /= 60.0
            return rt
    return 0.0


def _parse_spectrum(elem: ET.Element) -> Optional[MS1Spectrum | tuple]:
    """
    Parse a single <spectrum> element.

    Returns:
        MS1Spectrum             if this is an MS1 mass spectrum
        ('uv', rt, wl, int)     if this is an electromagnetic radiation spectrum
        None                    if it is neither
    """
    # Determine spectrum type
    is_ms1 = _has_cv(elem, _CV_MS1)
    is_uv  = _has_cv(elem, _CV_EM_RAD)

    if not is_ms1 and not is_uv:
        return None

    rt = _extract_rt(elem)

    # Decode binary arrays
    mz_arr = wl_arr = int_arr = None
    bda_list = _find(elem, "binaryDataArrayList")
    if bda_list is None:
        return None

    for bda in _iter(bda_list, "binaryDataArray"):
        if _has_cv(bda, _CV_MZ_ARRAY):
            mz_arr = _decode_binary_array(bda)
        elif _has_cv(bda, _CV_WAVELENGTH):
            wl_arr = _decode_binary_array(bda)
        elif _has_cv(bda, _CV_INTENSITY):
            int_arr = _decode_binary_array(bda)

    if is_ms1:
        if mz_arr is None or int_arr is None:
            return None
        polarity = "positive" if _has_cv(elem, _CV_POSITIVE) else "negative"
        return MS1Spectrum(rt=rt, polarity=polarity, mz=mz_arr, intensity=int_arr)

    if is_uv:
        if wl_arr is None or int_arr is None:
            return None
        return ("uv", rt, wl_arr, int_arr)

    return None


_WL_PATTERNS = [
    re.compile(r"Sig=(\d+)", re.IGNORECASE),
    re.compile(r"(\d+)\s*nm", re.IGNORECASE),
    re.compile(r"WVL[=:](\d+)", re.IGNORECASE),
    re.compile(r"WAVELENGTH[=:](\d+)", re.IGNORECASE),
]


def _wavelength_from_id(native_id: str) -> Optional[float]:
    for pat in _WL_PATTERNS:
        m = pat.search(native_id)
        if m:
            val = float(m.group(1))
            if 190 <= val <= 800:
                return val
    return None


def _parse_chromatogram(elem: ET.Element) -> Optional[UVChromatogram]:
    """
    Parse a <chromatogram> element.

    Returns a UVChromatogram only for DAD/UV chromatograms (identified by
    native ID keywords). Returns None for TIC and other non-UV traces.
    """
    native_id = elem.get("id", "")
    upper = native_id.upper()

    # Accept DAD, PDA, UV, WAVELENGTH in the ID; skip TIC and SRM traces
    if not any(kw in upper for kw in ("DAD", "PDA", "UV", "WAVELENGTH")):
        return None

    wavelength = _wavelength_from_id(native_id)

    bda_list = _find(elem, "binaryDataArrayList")
    if bda_list is None:
        return None

    time_arr = int_arr = None
    for bda in _iter(bda_list, "binaryDataArray"):
        # Time array: MS:1000595 (time array) or first non-intensity array
        if _has_cv(bda, "MS:1000595") or (time_arr is None and not _has_cv(bda, _CV_INTENSITY)):
            time_arr = _decode_binary_array(bda)
        elif _has_cv(bda, _CV_INTENSITY):
            int_arr = _decode_binary_array(bda)

    if time_arr is None or int_arr is None or len(time_arr) == 0:
        return None

    # Convert seconds → minutes if needed (chromatogram times are usually in seconds)
    # Heuristic: if max time > 60 min it's probably in seconds
    if time_arr.max() > 60:
        time_arr = time_arr / 60.0

    return UVChromatogram(
        native_id=native_id,
        wavelength=wavelength,
        times=time_arr,
        intensities=int_arr,
    )


# ── Public API ────────────────────────────────────────────────────────────────

def parse(filepath: str) -> ExperimentData:
    """
    Parse an mzML file and return all three data channels.

    Args:
        filepath: Path to the .mzML file.

    Returns:
        ExperimentData with ms1_spectra, uv_spectra, and uv_chromatograms populated.
        elsd_chromatogram is always None (mzML files never contain ELSD data).
    """
    result = ExperimentData(source_path=filepath, source_kind="mzml")
    uv_rows: list[tuple] = []  # (rt, wavelengths, intensities)

    # Use iterparse so large files do not require holding the full DOM
    context = ET.iterparse(filepath, events=("end",))
    for event, elem in context:
        # Strip namespace to get the local tag
        local = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag

        if local == "spectrum":
            parsed = _parse_spectrum(elem)
            if isinstance(parsed, MS1Spectrum):
                result.ms1_spectra.append(parsed)
            elif isinstance(parsed, tuple) and parsed[0] == "uv":
                _, rt, wl, intensity = parsed
                uv_rows.append((rt, wl, intensity))
            elem.clear()

        elif local == "chromatogram":
            chrom = _parse_chromatogram(elem)
            if chrom is not None:
                result.uv_chromatograms.append(chrom)
            elem.clear()

    # Assemble 2D DAD matrix if we collected UV scans
    if uv_rows:
        uv_rows.sort(key=lambda r: r[0])  # sort by RT
        times = np.array([r[0] for r in uv_rows])
        # Use wavelength array from first scan (assumed consistent across file)
        wavelengths = uv_rows[0][1]
        data = np.vstack([r[2] for r in uv_rows])
        result.uv_spectra = UVSpectraMatrix(
            times=times,
            wavelengths=wavelengths,
            data=data,
        )

    return result
