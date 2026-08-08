"""
Shared spectra data classes for MSCheck.

These are the container types returned by any input backend (mzml_parser,
rainbow_parser, etc.). Keeping them in a single module lets backends stay
independent and lets downstream modules (`mzspectrum`, `uvspectrum`,
`elsdchromatogram`) consume one canonical shape.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class MS1Spectrum:
    rt: float            # retention time in minutes
    polarity: str        # 'positive', 'negative', or 'unknown'
    mz: np.ndarray
    intensity: np.ndarray


@dataclass
class UVChromatogram:
    """Single-wavelength UV/DAD trace."""
    native_id: str
    wavelength: Optional[float]  # nm; None if not known
    times: np.ndarray            # minutes
    intensities: np.ndarray      # mAU or detector counts


@dataclass
class UVSpectraMatrix:
    """
    Full 2D DAD dataset: one UV spectrum per scan.

    Shape:
        times       : (n_scans,)                  minutes
        wavelengths : (n_wavelengths,)             nm
        data        : (n_scans, n_wavelengths)     mAU or detector counts
    """
    times: np.ndarray
    wavelengths: np.ndarray
    data: np.ndarray

    def at_wavelength(self, nm: float) -> tuple[np.ndarray, np.ndarray]:
        """Return (times, intensities) for the wavelength closest to *nm*."""
        idx = int(np.argmin(np.abs(self.wavelengths - nm)))
        return self.times, self.data[:, idx]

    def best_wavelength(self) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Return (wavelength_nm, times, intensities) for the wavelength with
        the largest integrated absorbance across the run.
        """
        areas = np.trapz(np.maximum(self.data, 0), self.times, axis=0)
        idx = int(np.argmax(areas))
        return float(self.wavelengths[idx]), self.times, self.data[:, idx]


@dataclass
class ELSDChromatogram:
    """
    Evaporative Light Scattering Detector trace.

    ELSD is a single-channel detector — one intensity value per scan.
    """
    times: np.ndarray            # minutes
    intensities: np.ndarray      # detector counts / mV (vendor-dependent)


@dataclass
class ExperimentData:
    """
    Parsed contents of one experiment (one mzML file or one vendor directory).

    Not all channels are always populated:
      • mzML files provide MS + UV (as available); never ELSD.
      • Agilent .D / Waters .raw via Rainbow provide MS + UV + ELSD (as available).
    """
    ms1_spectra: list[MS1Spectrum] = field(default_factory=list)
    uv_spectra: Optional[UVSpectraMatrix] = None
    uv_chromatograms: list[UVChromatogram] = field(default_factory=list)
    elsd_chromatogram: Optional[ELSDChromatogram] = None
    source_path: str = ""
    source_kind: str = ""  # 'mzml' | 'agilent' | 'waters'
