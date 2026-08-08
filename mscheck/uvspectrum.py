"""UV chromatogram and 2D DAD spectrum from mzML files or vendor directories."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from .spectra_types import ExperimentData, UVChromatogram, UVSpectraMatrix
from .reader import read as _read

logger = logging.getLogger(__name__)


class UVSpectrum:
    """
    UV data from one mzML file or one Agilent/Waters vendor directory.

    Prioritises the 2D DAD spectral matrix (electromagnetic radiation spectra,
    typically 190–400 nm at every scan) when present.  Falls back to any
    single-wavelength chromatograms.

    Public attributes
    -----------------
    has_uv : bool
    uv_spectra : UVSpectraMatrix or None   — 2D DAD matrix (n_scans × n_wavelengths)
    uv_chromatograms : list[UVChromatogram] — single-wavelength traces (fallback)
    """

    def __init__(self, filepath: str) -> None:
        self._filepath = filepath
        self._data: ExperimentData = _read(filepath)

        self.uv_spectra: Optional[UVSpectraMatrix] = self._data.uv_spectra
        self.uv_chromatograms: list[UVChromatogram] = self._data.uv_chromatograms
        self.has_uv = (self.uv_spectra is not None) or bool(self.uv_chromatograms)

        if self.uv_spectra is not None:
            logger.info(
                "Loaded 2D DAD matrix: %d scans x %d wavelengths (%.0f-%.0f nm)",
                self.uv_spectra.data.shape[0],
                self.uv_spectra.data.shape[1],
                self.uv_spectra.wavelengths[0],
                self.uv_spectra.wavelengths[-1],
            )
        elif self.uv_chromatograms:
            wls = [c.wavelength for c in self.uv_chromatograms if c.wavelength]
            logger.info("Loaded %d UV chromatogram(s): %s nm",
                        len(self.uv_chromatograms), wls)
        else:
            logger.info("No UV data found in %s", filepath)

    # ── Primary access ────────────────────────────────────────────────────────

    def get_chromatogram(self, target_wavelength: Optional[float] = None
                         ) -> tuple[np.ndarray, np.ndarray]:
        """
        Return (times_min, intensities) for *target_wavelength* nm.

        Uses the 2D DAD matrix when available.  Falls back to single-wavelength
        chromatograms otherwise.  If *target_wavelength* is None, returns the
        chromatogram with the highest integrated area across the run.
        """
        if not self.has_uv:
            return np.array([]), np.array([])

        if self.uv_spectra is not None:
            if target_wavelength is None:
                _, times, intensities = self.uv_spectra.best_wavelength()
            else:
                times, intensities = self.uv_spectra.at_wavelength(target_wavelength)
            return times, intensities

        chrom = self._closest_chromatogram(target_wavelength)
        if chrom is None:
            return np.array([]), np.array([])
        return chrom.times, chrom.intensities

    def get_best_wavelength(self) -> Optional[float]:
        """Return the wavelength (nm) with the largest integrated peak area."""
        if self.uv_spectra is not None:
            wl, _, _ = self.uv_spectra.best_wavelength()
            return wl
        if self.uv_chromatograms:
            best = max(self.uv_chromatograms, key=lambda c: c.intensities.max())
            return best.wavelength
        return None

    def get_wavelengths(self) -> list[float]:
        """Return all available wavelengths (nm)."""
        if self.uv_spectra is not None:
            return self.uv_spectra.wavelengths.tolist()
        return [c.wavelength for c in self.uv_chromatograms if c.wavelength is not None]

    def get_multi_wavelength_data(self) -> dict:
        """
        Return UV data for every wavelength as a dict of parallel lists.

        Keys:
            'wavelengths'     : list[float] (nm)
            'retention_times' : list[np.ndarray] (minutes), one per wavelength
            'intensities'     : list[np.ndarray], one per wavelength
        """
        if self.uv_spectra is not None:
            n = len(self.uv_spectra.wavelengths)
            return {
                "wavelengths": self.uv_spectra.wavelengths.tolist(),
                "retention_times": [self.uv_spectra.times] * n,
                "intensities": [self.uv_spectra.data[:, i] for i in range(n)],
            }

        return {
            "wavelengths": [c.wavelength for c in self.uv_chromatograms],
            "retention_times": [c.times for c in self.uv_chromatograms],
            "intensities": [c.intensities for c in self.uv_chromatograms],
        }

    def get_max_absorption_chromatogram(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Return (times_min, intensities) where each point is the maximum
        absorbance across all wavelengths at that scan time.
        """
        if not self.has_uv:
            return np.array([]), np.array([])

        if self.uv_spectra is not None:
            return self.uv_spectra.times, np.max(self.uv_spectra.data, axis=1)

        if not self.uv_chromatograms:
            return np.array([]), np.array([])
        times = self.uv_chromatograms[0].times
        stacked = np.vstack([c.intensities for c in self.uv_chromatograms])
        return times, np.max(stacked, axis=0)

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _closest_chromatogram(self, wavelength: Optional[float]) -> Optional[UVChromatogram]:
        if not self.uv_chromatograms:
            return None
        if wavelength is None:
            return self.uv_chromatograms[0]
        known = [(c, c.wavelength) for c in self.uv_chromatograms if c.wavelength is not None]
        if not known:
            return self.uv_chromatograms[0]
        return min(known, key=lambda x: abs(x[1] - wavelength))[0]
