"""Mass spectrum from .mxML files"""
from pyopenms import *


class MassSpectrum(object):
    """
    Creates mass specturm object using the pyopenms package
    """

    def __init__(self, mzMLfilepath: str, mode: str = "Positive"):
        """
        MassSpecturm constructor
        Args:
            mzMLfilepath (str): path to .mxML file
            mode (str): specify polarity ("Positive" or "Negative") of mass spctra mode to search ions for. Default set
                        to Positive.
        """
        self._filepath = mzMLfilepath
        self._mode = mode
        self._exp = self._create_openms()
        self.MSdata = self._get_ms_data()

    def _create_openms(self) -> None:
        """
        Creates an openms class object
        Returns:
            openmsobject: openms object created using .mxML filepath
        """
        openms = MSExperiment()
        MzMLFile().load(self._filepath, openms)
        if self._mode == "Positive":
            openms = [spectrum for i, spectrum in enumerate(openms) if i % 2 == 0]
        if self._mode == "Negative":
            openms = [spectrum for i, spectrum in enumerate(openms) if i % 2 != 0]
        return openms

    def _get_ms_data(self) -> dict:
        """
        Creates a dictionary of MS data
        Returns:
            MS_data (dict): a dictionary of mass spectrum peaks, peak instensity, total ion count
                            and retention time data for a spectrum
        """
        mz_data = [spectrum.get_peaks() for spectrum in self._exp]
        TIC_data = [round(sum(data[1])) for data in mz_data]
        RT_data = [spectrum.getRT() for spectrum in self._exp]

        MS_data = {
            "mz_data": mz_data,
            "TIC": np.array(TIC_data),
            "RT": np.array(RT_data) / 60,
        }

        return MS_data
