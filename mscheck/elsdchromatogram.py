"""ELSD chromatogram from vendor binary directories.

ELSD (Evaporative Light Scattering Detector) is a single-channel detector.
Its data can only be read from vendor binary files (Agilent .D, Waters .raw)
via the Rainbow parser — msConvert does not export ELSD channels into mzML.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

from .spectra_types import ELSDChromatogram, ExperimentData
from .reader import read as _read

logger = logging.getLogger(__name__)


class ELSDChromatogram:
    """
    ELSD trace from one Agilent/Waters vendor directory.

    Public attributes
    -----------------
    has_elsd : bool
    trace : ELSDChromatogram or None
    """

    def __init__(self, filepath: str) -> None:
        self._filepath = filepath
        self._data: ExperimentData = _read(filepath)
        self.trace: Optional[ELSDChromatogram] = self._data.elsd_chromatogram
        self.has_elsd = self.trace is not None

        if self.has_elsd:
            assert self.trace is not None
            logger.info(
                "Loaded ELSD trace: %d points, RT %.3f–%.3f min",
                len(self.trace.times), self.trace.times.min(), self.trace.times.max(),
            )
        else:
            logger.info("No ELSD data in %s", filepath)

    # ── Access ────────────────────────────────────────────────────────────────

    def get_chromatogram(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (times_min, intensities). Empty arrays if no ELSD data."""
        if self.trace is None:
            return np.array([]), np.array([])
        return self.trace.times, self.trace.intensities

    def peak_area(self, rt_start: float, rt_end: float) -> float:
        """
        Integrate the ELSD signal between *rt_start* and *rt_end* (minutes).

        Returns 0.0 if no ELSD data or the window contains no points.
        """
        if self.trace is None:
            return 0.0
        t, y = self.trace.times, self.trace.intensities
        mask = (t >= rt_start) & (t <= rt_end)
        if not mask.any():
            return 0.0
        return float(np.trapz(y[mask], t[mask]))
