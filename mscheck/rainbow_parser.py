"""
Rainbow-based reader for Agilent (.D, .dx) and Waters (.raw) vendor directories.

Uses the ``rainbow-api`` package to read binary vendor files directly, providing:
  • ELSD chromatograms (unavailable via mzML)
  • 2D DAD spectral matrix (as with mzML, but read from the native binary)
  • MS1 spectra

Rainbow is an optional dependency. If it is not installed, importing this module
raises ImportError with a helpful message. Callers should catch this and fall back
to the mzML backend.

Polarity limitation
-------------------
Rainbow's Agilent/Waters MS parsers do not expose per-scan polarity information.
When the input directory contains a single-polarity acquisition, this is fine.
For mixed-polarity (fast-switching) acquisitions, all MS1 spectra returned from
this parser are labelled ``polarity='unknown'`` and the caller (``MassSpectrum``)
treats them as matching whichever mode was requested.
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np

from .spectra_types import (
    ELSDChromatogram,
    ExperimentData,
    MS1Spectrum,
    UVChromatogram,
    UVSpectraMatrix,
)

try:
    import rainbow as _rb  # type: ignore
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "rainbow_parser requires the 'rainbow-api' package. "
        "Install with: pip install rainbow-api"
    ) from exc


_VENDOR_EXTENSIONS = {
    ".d": "agilent",
    ".dx": "agilent",
    ".raw": "waters",
}


def is_vendor_directory(path: str) -> bool:
    """Return True if *path* is a directory that Rainbow can read."""
    if not os.path.isdir(path):
        return False
    ext = os.path.splitext(path)[1].lower()
    return ext in _VENDOR_EXTENSIONS


def _detect_kind(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    return _VENDOR_EXTENSIONS.get(ext, "unknown")


# ── DataFile → MSCheck type converters ────────────────────────────────────────

def _ms_datafile_to_spectra(df, requested_polarity: Optional[str] = None
                            ) -> list[MS1Spectrum]:
    """
    Convert a rainbow MS DataFile (2D matrix of RT × m/z intensities) into a
    list of MS1Spectrum. Zero-intensity m/z channels are dropped per scan.

    Rainbow does not expose per-scan polarity; all spectra receive
    polarity ``requested_polarity`` if given, else ``'unknown'``.
    """
    times = np.asarray(df.xlabels, dtype=np.float64)  # already minutes
    mz_labels = np.asarray(df.ylabels, dtype=np.float64)
    data = np.asarray(df.data)  # shape (n_scans, n_mz)

    polarity = requested_polarity or "unknown"

    spectra: list[MS1Spectrum] = []
    for i in range(data.shape[0]):
        row = data[i]
        nz = np.flatnonzero(row)
        spectra.append(MS1Spectrum(
            rt=float(times[i]),
            polarity=polarity,
            mz=mz_labels[nz].copy(),
            intensity=row[nz].astype(np.float64, copy=True),
        ))
    return spectra


def _uv_datafile_to_matrix(df) -> Optional[UVSpectraMatrix]:
    """Convert a rainbow UV DataFile (RT × wavelength) into a UVSpectraMatrix."""
    times = np.asarray(df.xlabels, dtype=np.float64)
    wavelengths = np.asarray(df.ylabels, dtype=np.float64)
    data = np.asarray(df.data, dtype=np.float64)

    if data.size == 0 or wavelengths.size == 0:
        return None
    return UVSpectraMatrix(times=times, wavelengths=wavelengths, data=data)


def _uv_datafile_to_chromatogram(df) -> Optional[UVChromatogram]:
    """Convert a single-wavelength UV DataFile (ylabels length 1) into a UVChromatogram."""
    times = np.asarray(df.xlabels, dtype=np.float64)
    data = np.asarray(df.data, dtype=np.float64)
    wavelengths = np.asarray(df.ylabels, dtype=np.float64)

    if data.ndim == 2 and data.shape[1] == 1:
        intensity = data[:, 0]
    else:
        intensity = data.ravel()

    if intensity.size == 0:
        return None

    wavelength: Optional[float] = float(wavelengths[0]) if wavelengths.size else None
    return UVChromatogram(
        native_id=df.name,
        wavelength=wavelength,
        times=times,
        intensities=intensity,
    )


def _elsd_datafile_to_chromatogram(df) -> Optional[ELSDChromatogram]:
    """Convert a rainbow ELSD DataFile (single channel) into an ELSDChromatogram."""
    times = np.asarray(df.xlabels, dtype=np.float64)
    data = np.asarray(df.data, dtype=np.float64)

    if data.ndim == 2 and data.shape[1] == 1:
        intensity = data[:, 0]
    else:
        intensity = data.ravel()

    if intensity.size == 0 or times.size == 0:
        return None
    return ELSDChromatogram(times=times, intensities=intensity)


# ── Public API ────────────────────────────────────────────────────────────────

def parse(path: str, requested_polarity: Optional[str] = None) -> ExperimentData:
    """
    Read an Agilent (.D, .dx) or Waters (.raw) directory via Rainbow.

    Args:
        path: Directory path.
        requested_polarity: If given ('positive' or 'negative'), all MS1 spectra
            are tagged with this polarity. Rainbow does not currently expose
            per-scan polarity for Agilent/Waters; use only when you know the
            acquisition was single-polarity.

    Returns:
        ExperimentData with whichever channels the vendor file contained.
    """
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Not a directory: {path}")

    dd = _rb.read(path)

    result = ExperimentData(source_path=path, source_kind=_detect_kind(path))

    ms_files   = dd.by_detector.get("MS",   [])
    uv_files   = dd.by_detector.get("UV",   [])
    elsd_files = dd.by_detector.get("ELSD", [])

    # ── MS: pick the first MS DataFile ──
    if ms_files:
        result.ms1_spectra = _ms_datafile_to_spectra(ms_files[0], requested_polarity)

    # ── UV: prefer a 2D matrix; also collect any single-wavelength traces ──
    for uf in uv_files:
        n_wavelengths = int(np.asarray(uf.ylabels).size)
        if n_wavelengths > 1 and result.uv_spectra is None:
            result.uv_spectra = _uv_datafile_to_matrix(uf)
        else:
            chrom = _uv_datafile_to_chromatogram(uf)
            if chrom is not None:
                result.uv_chromatograms.append(chrom)

    # ── ELSD: pick the first ELSD DataFile ──
    if elsd_files:
        result.elsd_chromatogram = _elsd_datafile_to_chromatogram(elsd_files[0])

    return result
