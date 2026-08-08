"""Mass spectrum from mzML files or vendor directories."""

from __future__ import annotations

import numpy as np

from .reader import read as _read


class MassSpectrum:
    """
    Loads MS1 spectra from an mzML file or vendor directory, filtered by polarity.

    Public attributes
    -----------------
    mode : str        — "Positive" or "Negative"
    MSdata : dict     — {'mz_data': list of (mz, intensity) tuples,
                         'TIC':     np.ndarray of rounded per-scan sums,
                         'RT':      np.ndarray of retention times in minutes}
    """

    def __init__(self, filepath: str, mode: str = "Positive") -> None:
        self._filepath = filepath
        self.mode = mode
        self.MSdata = self._load(filepath, mode)

    @staticmethod
    def _load(path: str, mode: str) -> dict:
        target = "positive" if mode == "Positive" else "negative"
        data = _read(path, requested_polarity=target)

        # 'unknown' polarity comes from Rainbow-loaded vendor data — accept it
        # under either mode. mzML data always has an explicit polarity.
        spectra = [
            s for s in data.ms1_spectra
            if s.polarity == target or s.polarity == "unknown"
        ]

        mz_data = [(s.mz, s.intensity) for s in spectra]
        tic = np.array([round(float(s.intensity.sum())) for s in spectra])
        rt = np.array([s.rt for s in spectra])

        return {"mz_data": mz_data, "TIC": tic, "RT": rt}
