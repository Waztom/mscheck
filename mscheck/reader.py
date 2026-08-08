"""
Input dispatcher for MSCheck.

Given a path, selects the appropriate parser:
    *.mzML                  → mzml_parser (built-in, stdlib only)
    *.D / *.dx / *.raw (dir)→ rainbow_parser (optional dep: rainbow-api)

Downstream modules (`MassSpectrum`, `UVSpectrum`, `ELSDChromatogram`) call
``read()`` and are agnostic to the underlying format.
"""

from __future__ import annotations

import os
from typing import Optional

from .spectra_types import ExperimentData
from . import mzml_parser
from . import rainbow_parser

_RAINBOW_EXTENSIONS = {".d", ".dx", ".raw"}


def _looks_like_vendor_directory(path: str) -> bool:
    if not os.path.isdir(path):
        return False
    return os.path.splitext(path)[1].lower() in _RAINBOW_EXTENSIONS


def read(path: str, requested_polarity: Optional[str] = None) -> ExperimentData:
    """
    Read an mzML file or a vendor directory into a common ``ExperimentData``.

    Args:
        path: Path to a .mzML file, or to an Agilent .D / .dx / Waters .raw
            directory.
        requested_polarity: Only relevant when reading a vendor directory via
            Rainbow (which does not expose per-scan polarity). Ignored for mzML.

    Raises:
        ImportError: if a vendor directory is passed but ``rainbow-api`` is not
            installed.
        FileNotFoundError: if the path does not exist.
    """
    if _looks_like_vendor_directory(path):
        return rainbow_parser.parse(path, requested_polarity=requested_polarity)

    if not os.path.isfile(path):
        raise FileNotFoundError(f"Not a file or recognised vendor directory: {path}")

    return mzml_parser.parse(path)
