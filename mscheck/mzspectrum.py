"""Mass spectrum from .mzML files."""

from __future__ import annotations

import numpy as np

from .mzml_parser import parse as _parse_mzml


class MassSpectrum:
    """
    Loads MS1 spectra from an mzML file, filtered by polarity.

    Public attributes
    -----------------
    mode : str        — "Positive" or "Negative"
    MSdata : dict     — {'mz_data': list of (mz, intensity) tuples,
                         'TIC':     np.ndarray of rounded per-scan sums,
                         'RT':      np.ndarray of retention times in minutes}
    """

    def __init__(self, mzMLfilepath: str, mode: str = "Positive") -> None:
        self._filepath = mzMLfilepath
        self.mode = mode
        self.MSdata = self._load(mzMLfilepath, mode)

    @staticmethod
    def _load(mzMLfilepath: str, mode: str) -> dict:
        data = _parse_mzml(mzMLfilepath)

        target = "positive" if mode == "Positive" else "negative"
        spectra = [s for s in data.ms1_spectra if s.polarity == target]

        mz_data = [(s.mz, s.intensity) for s in spectra]
        tic = np.array([round(float(s.intensity.sum())) for s in spectra])
        rt = np.array([s.rt for s in spectra])

        return {"mz_data": mz_data, "TIC": tic, "RT": rt}
