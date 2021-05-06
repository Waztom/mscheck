"""Analyse spectrum class """
from __future__ import annotations
from scipy.signal import find_peaks, peak_widths
import numpy as np
from .utils import get_smiles, get_mol, get_MW, get_path_leaf
from .report import create_report_plot
from .spectrum import MassSpectrum


class AnalyseSpectrum(MassSpectrum):
    """
    Analyses the MassSpectrum class object
    """

    def __init__(self, mzMLfilepath: str, mode: str = "Positive"):
        """
        Analyse spectrum constructor
        mzMLfilepath (str): path to .mxML file
        mode (str): specify polarity ("Positive" or "Negative") of mass spctra mode to search ions for. Default set
                    to Positive.
        """
        super().__init__(mzMLfilepath, mode)
        self.MSpeakdata = self._get_ms_peak_data()

    def _get_ms_peak_data(self) -> dict:
        """
        Uses scipy peak finding and peak width algorithims to find peak indexes and widths
        to gather MS data from
        Returns:
            MS_peak_data dict: a dictionary of mass to charge ratio (mz), total ion count and the
                            retention time at points where peaks have been found
        """
        MS_peak_RT = []
        MS_peak_TIC = []
        MS_peak_mz_max = []
        MS_peak_mz_data = []

        peak_indices, peak_properties = find_peaks(self.MSdata["TIC"])
        peak_width_results = peak_widths(self.MSdata["TIC"], peak_indices, rel_height=0.5)

        if len(peak_indices) > 0:
            for peak_index, FWHM_height in zip(
                peak_indices,
                peak_width_results[1],
            ):
                FWHM_indices = self.get_FWHM_indices(FWHM_height, peak_index)
                MS_peak_RT.append(
                    [self.MSdata["RT"][FWHM_index] for FWHM_index in FWHM_indices]
                )
                MS_peak_TIC.append(
                    [self.MSdata["TIC"][FWHM_index] for FWHM_index in FWHM_indices]
                )

                peak_mz_data = [
                    self.MSdata["mz_data"][FWHM_index] for FWHM_index in FWHM_indices
                ]

                MS_peak_mz_data.append(peak_mz_data)

                MS_peak_mz_max.append([self.get_max_mz(data) for data in peak_mz_data])

        return {
            "RT": np.array([item for sublist in MS_peak_RT for item in sublist]),
            "TIC": np.array([item for sublist in MS_peak_TIC for item in sublist]),
            "mz_data": [item for sublist in MS_peak_mz_data for item in sublist],
            "mz_max": np.array([item for sublist in MS_peak_mz_max for item in sublist]),
        }

    def get_FWHM_indices(self, FWHM_height: float, peak_index: int) -> list:
        """
        Searches for indices left and right of the peak index for total ion count
        values above the height from the peaks full width half maximum (FWHM).
        This enables searching of impure signals with mixed compounds.
        Args:
            FWHM_height (float): height of the FWHM threshold
            peak_index (int): index of the peak found

        """
        indices = []
        max_TIC = self.MSdata["TIC"][peak_index]

        left_indices = []
        left_index = peak_index
        left_TIC = self.MSdata["TIC"][peak_index - 1]

        while left_TIC > FWHM_height and left_TIC < max_TIC:
            left_index = left_index - 1
            if left_index >= 0:
                left_TIC = self.MSdata["TIC"][left_index]
                left_indices.append(left_index)

        indices.extend(left_indices)
        indices.append(peak_index)

        right_indices = []
        right_index = peak_index
        right_TIC = self.MSdata["TIC"][peak_index + 1]

        while right_TIC > FWHM_height and right_TIC < max_TIC:
            right_index = right_index + 1
            if right_index <= len(self.MSdata["TIC"]):
                right_TIC = self.MSdata["TIC"][right_index]
                right_indices.append(right_index)

        indices.extend(right_indices)

        return sorted(indices)

    def get_max_mz(self, mz_data: list) -> int or None:
        """
        finds the maximum mz value from a list
        Args:
            mz_data (list): list of two list containing mz
                            and intensity values from fragment pattern
        """

        mz_values = mz_data[0]
        intensity_values = mz_data[1]
        if intensity_values.size:
            mz_max = int(round(mz_values[np.argmax(intensity_values)]))
            return mz_max
        else:
            return None

    def get_match_indices(self, mass_ion: int) -> list:
        """
        Finds indices where mass of ions found in mz data
        Args:
            mass_ion (int): rounded molecular weight of the target
                            compound + ionisation ion
        """
        match_test = list(self.MSpeakdata["mz_max"]).count(mass_ion)

        if match_test != 0:
            match_indices = [
                i for i, x in enumerate(self.MSpeakdata["mz_max"]) if x == mass_ion
            ]
            return mass_ion, match_indices
        else:
            return mass_ion, []

    def analyse(self, compoundsmiles: str, ionstoadd: list, tolerance: int) -> dict:
        """
        Performs analysis of spectrum
        Args:
            compoundsmiles (str): SMILES string for the target compound ebing analysed
            ionstoadd (list): list of ions as neutral SMILES to add to target compound mass for
                              searching mass spectrum eg. ["[H]", "[Na]"]
            tolerance (int): tolerance set for finding a match. Eg. tolerance set to 1 will
                             search for target compound with mass in spectrum of:
                             MW target compound plus/minus 1
        """
        self.compound_mol = get_mol(compoundsmiles)
        self.compound_MW = get_MW(self.compound_mol)
        ions_matched = []
        RT_matched = []
        TIC_matched = []
        mz_data_matched = []
        mz_strongest_value = []

        ions_to_add_mols = [get_mol(ion) for ion in ionstoadd]
        ions_to_add_MW = [get_MW(mol) for mol in ions_to_add_mols]

        for ion_to_add_mol, MW in zip(ions_to_add_mols, ions_to_add_MW):
            parent_mass = MW + self.compound_MW
            test_ion_masses = [
                parent_mass - tolerance,
                parent_mass,
                parent_mass + tolerance,
            ]

            match_results = [self.get_match_indices(ion) for ion in test_ion_masses]

            for result in match_results:
                if result[1]:
                    ions_matched.append((get_smiles(ion_to_add_mol), result[0]))
                    RT_matched.append(self.MSpeakdata["RT"][result[1]])
                    TIC_matched.append(self.MSpeakdata["TIC"][result[1]])
                    mz_data_matched.append(
                        [self.MSpeakdata["mz_data"][index] for index in result[1]]
                    )
                    mz_strongest_value.append(
                        self.get_strongest_mz_pattern(
                            [self.MSpeakdata["mz_data"][index] for index in result[1]]
                        )
                    )

        self.Matchdata = {
            "ions": ions_matched,
            "RT": RT_matched,
            "TIC": TIC_matched,
            "mz_data": mz_data_matched,
            "mz_strongest": mz_strongest_value,
        }

    def get_strongest_mz_pattern(self, mz_values: list) -> largest_mz_pattern:
        """
        Finds largest mz signal from list of mz patterns and returns mz maz,
        the intensity and index of max value
            Args:
                mz_data (list): lists of list with two lists of mz and intensity values
                                for an ion fragment pattern
        """
        mz_intensities = [mz_value[1] for mz_value in mz_values]
        mz_masses = [mz_value[0] for mz_value in mz_values]

        for index, intensity_values in enumerate(mz_intensities):
            if index == 0:
                max_value = np.amax(intensity_values)
                max_index = index
            else:
                max_value_test = np.amax(intensity_values)
                if max_value_test > max_value:
                    max_value = max_value_test
                    max_index = index
                else:
                    pass
        mz_masses_max = mz_masses[max_index]
        mz_intensities_max = mz_intensities[max_index]
        return mz_masses_max, mz_intensities_max, max_index

    def create_report(self, compound_name: str = None) -> MassCheckReport:
        no_plots = len(self.Matchdata["ions"])
        if not compound_name:
            compound_name = get_path_leaf(self._filepath)

        create_report_plot(
            RT_values=self.MSdata["RT"],
            TIC_values=self.MSdata["TIC"],
            no_plots=no_plots,
            mol=self.compound_mol,
            match_data=self.Matchdata,
            compound_name=compound_name,
        )
