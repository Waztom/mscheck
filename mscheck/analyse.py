#%%
"""Analyse spectrum class """
from scipy.signal import find_peaks, peak_widths
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, AnnotationBbox, OffsetImage
from matplotlib.backends.backend_pgf import PdfPages

from utils import get_smiles, get_mol, get_MW, get_molecule_image
from spectrum import MassSpectrum


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
        MS_peak_mz = []

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
                MS_peak_values = [
                    self._exp[FWHM_index].get_peaks() for FWHM_index in FWHM_indices
                ]

                MS_peak_mz.append([self.get_max_mz(value) for value in MS_peak_values])

        return {
            "RT": np.array([item for sublist in MS_peak_RT for item in sublist]),
            "TIC": np.array([item for sublist in MS_peak_TIC for item in sublist]),
            "mz": np.array([item for sublist in MS_peak_mz for item in sublist]),
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

    def get_max_mz(self, MS_values: list) -> int or None:
        """
        finds the maximum mz value from a list
        Args:
            MS_values (list): list of mz values from fragment pattern
        """

        mz_values = MS_values[0]
        intensity_values = MS_values[1]

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
        match_test = list(self.MSpeakdata["mz"]).count(mass_ion)
        if match_test != 0:
            match_indices = [
                i for i, x in enumerate(self.MSpeakdata["mz"]) if x == mass_ion
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

        self.match_data = {
            "ions_matched": ions_matched,
            "RT_matched": RT_matched,
            "TIC_matched": TIC_matched,
        }

        return self.match_data

    def get_report_figure(self):
        fig, ax = plt.subplots(figsize=(10.694, 7.264))
        pdf_pages = PdfPages("report.pdf")
        ax.plot(self.MSdata["RT"], self.MSdata["TIC"], color="k", zorder=-1)

        for ion_found, RT_values, TIC_values in zip(
            self.match_data["ions_matched"],
            self.match_data["RT_matched"],
            self.match_data["TIC_matched"],
        ):
            ax.scatter(RT_values, TIC_values, s=10.0, marker="x", color="r", zorder=1)
            TIC_max_index = np.argmax(TIC_values)
            xy = (RT_values[TIC_max_index], TIC_values[TIC_max_index])

            offsetbox = TextArea(
                "Ion: {} \nM+: {}\nRT: {}".format(
                    ion_found[0], ion_found[1], np.round(RT_values[TIC_max_index], 2)
                )
            )

            ab = AnnotationBbox(
                offsetbox,
                xy,
                xybox=(20, 30),
                xycoords="data",
                boxcoords=("offset points"),
                box_alignment=(0.0, 0.5),
                arrowprops=dict(arrowstyle="->"),
            )
            ax.add_artist(ab)

        compound_image = get_molecule_image(self.compound_mol)

        imagebox = OffsetImage(compound_image, zoom=0.4)
        xy = (np.amax(self.MSdata["RT"]), np.amax(self.MSdata["TIC"]))
        ab = AnnotationBbox(
            imagebox,
            xy,
            pad=0.1,
        )
        imagebox.image.axes = ax

        ax.add_artist(ab)
        plt.xlabel("Retention time (min)")
        plt.ylabel("Total ion count (TIC)")
        plt.title(self._filepath.split("/")[-1])
        pdf_pages.savefig(fig, dpi=200)
        pdf_pages.close()
        return plt.show()


analysis_test = AnalyseSpectrum(
    "/home/warren/XChem_projects/carmc/carmc/1AB-1001.mzML", mode="Positive"
)

test = analysis_test.analyse(
    compoundsmiles="O=C(c1ccco1)N1CCN(C(=O)N2CCN(c3ccccc3)CC2)CC1",
    ionstoadd=["[H]", "[Na]"],
    tolerance=1,
)


analysis_test = AnalyseSpectrum(
    "/home/warren/XChem_projects/carmc/carmc/1AB-1001.mzML", mode="Positive"
)

test = analysis_test.analyse(
    compoundsmiles="O=C(c1ccco1)N1CCN(C(=O)N2CCN(c3ccccc3)CC2)CC1",
    ionstoadd=["[H]", "[Na]"],
    tolerance=1,
)

analysis_test.get_report_figure()
