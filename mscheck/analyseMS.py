"""Analyse spectrum class"""

from __future__ import annotations
from scipy.signal import find_peaks, peak_widths
import numpy as np
import logging
from .utils import get_smiles, get_mol, get_MW, get_path_leaf
from .report import MSReport
from .mzspectrum import MassSpectrum


class AnalyseMS(MassSpectrum):
    """
    Analyses the MassSpectrum class object
    """

    def __init__(self, filepath: str, mode: str = "Positive"):
        """
        Analyse spectrum constructor
        filepath (str): path to .mzML file or vendor directory (.D, .dx, .raw)
        mode (str): specify polarity ("Positive" or "Negative") of mass spctra mode to search ions for. Default set
                    to Positive.
        """
        super().__init__(filepath, mode)
        self.MSpeakdata = self._get_ms_peak_data()
        self.logger = logging.getLogger("AnalyseMS")

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
        peak_width_results = peak_widths(
            self.MSdata["TIC"], peak_indices, rel_height=0.5
        )

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
            "mz_max": np.array(
                [item for sublist in MS_peak_mz_max for item in sublist]
            ),
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

    def get_max_match_indices(self, mass_ion: int) -> list:
        """
        Finds data indices where mass of ions matches with maximum signal
        from  mz pattern data ie a veyr dominant signal
        Args:
            mass_ion (int): rounded molecular weight of the target
                            compound + ionisation ion
        """
        match_test = list(self.MSpeakdata["mz_max"]).count(mass_ion)

        if match_test != 0:
            match_indices = [
                i for i, x in enumerate(self.MSpeakdata["mz_max"]) if x == mass_ion
            ]
            max_mz_match = True
            return mass_ion, match_indices, max_mz_match
        else:
            max_mz_match = False
            return mass_ion, [], max_mz_match

    def get_any_match_indices(self, mass_ion: int) -> list:
        """
        Finds any data indices where mass of ions found in mz data
        ie a secondary or less dominant signal in the mz pattern
        Args:
            mass_ion (int): rounded molecular weight of the target
                            compound + ionisation ion
        """

        match_indices = []
        for i, mz_data in enumerate(self.MSpeakdata["mz_data"]):
            for mass in mz_data[0]:
                if round(mass) == mass_ion:
                    match_indices.append(i)
        if match_indices:
            max_mz_match = False
            return mass_ion, match_indices, max_mz_match
        else:
            max_mz_match = False
            return mass_ion, [], max_mz_match

    def get_extracted_ion_count(
        self, RT_data: list, mz_data: list, ion_mass: int
    ) -> list:
        """
        Extracts ion data from mz_data list
        Args:
            mz_data (list): list of mz data
            ion_mass (int): Masses ion to extract and match from mz_data
        Returns:
            extracted_ion_data (list): list of intensity values matched to ion to extract
        """
        ion_masses = [ion_mass, ion_mass + 1, ion_mass - 1]
        EIC_data = []
        for e, mz in enumerate(mz_data):
            extracted_ion_signal = 0
            for i, mass in enumerate(mz[0]):
                if round(mass) in ion_masses:
                    extracted_ion_signal += mz[1][i]
            EIC_data.append((RT_data[e], extracted_ion_signal))
        return EIC_data

    def analyse(
        self,
        compoundsmiles: str,
        ionstoadd: list,
        tolerance: int,
        ionstosub: list = None,
        custom_mw: float = None,
    ) -> dict:
        """
        Performs analysis of spectrum
        Args:
            compoundsmiles (str): SMILES string for the target compound being analysed
            ionstoadd (list): list of ions as neutral SMILES to add to target compound mass for
                              searching mass spectrum eg. ["[H]", "[Na]"]
            ionstosub (list): list of ions as neutral SMILES to remove from a target compound
                              mass for searching mass spectrum
                              eg. ["CC(C)=C","O=C=O.CC(C)=C"] for -tBu and -Boc fragmentation
            tolerance (int): tolerance set for finding a match. Eg. tolerance set to 1 will
                             search for target compound with mass in spectrum of:
                             MW target compound plus/minus 1
            custom_mw (float): Optional custom molecular weight to use for analysis. If provided,
                               this will override the molecular weight calculated from the SMILES
                               string.
        """
        self.logger.info(f"Analyzing compound with SMILES: {compoundsmiles}")
        self.logger.info(f"Ionization mode: {self.mode}, Tolerance: {tolerance}")
        self.logger.info(f"Ions to add: {ionstoadd}")
        if ionstosub:
            self.logger.info(f"Ions to subtract: {ionstosub}")

        # Parse SMILES to molecule
        self.compound_mol = get_mol(compoundsmiles)
        if self.compound_mol is None:
            self.logger.error(f"Failed to parse SMILES: {compoundsmiles}")
            return {}

        # Use custom molecular weight if provided, otherwise calculate from SMILES
        if custom_mw is not None:
            self.compound_MW = custom_mw
            self.logger.info(f"Using custom molecular weight: {custom_mw:.4f} Da")
        else:
            self.compound_MW = get_MW(self.compound_mol)
            self.logger.info(
                f"Using calculated molecular weight: {self.compound_MW:.4f} Da"
            )

        # Initialize data storage
        EIC_data = []
        max_mz_match = []
        ions_matched = []
        RT_matched = []
        TIC_matched = []
        mz_data_matched = []
        mz_strongest_value = []

        # Process ions to add
        self.logger.info(f"Processing {len(ionstoadd)} ions to add")
        ions_to_add_mols = [get_mol(ion) for ion in ionstoadd]
        ions_to_add_MW = [get_MW(mol) for mol in ions_to_add_mols]

        # Process ions to subtract (if any)
        if ionstosub:
            self.logger.info(f"Processing {len(ionstosub)} ions to subtract")
            ions_to_sub_mols = [get_mol(ion) for ion in ionstosub]
            # Mass +1 to account for resultant +H adduct observed - assuming positive ionisation method
            ions_to_sub_MW = [1 - get_MW(mol) for mol in ions_to_sub_mols]
        else:
            ions_to_sub_mols = []
            ions_to_sub_MW = []

        # Combine all ion modifications
        ions_to_alter_mols = ions_to_add_mols + ions_to_sub_mols
        ions_to_alter_MW = ions_to_add_MW + ions_to_sub_MW

        # Process each ion
        for ion_index, (ion_to_alter_mol, MW) in enumerate(
            zip(ions_to_alter_mols, ions_to_alter_MW)
        ):
            # Calculate expected m/z based on mode
            if self.mode == "Positive":
                parent_mass = MW + self.compound_MW
                self.logger.info(
                    f"Ion {ion_index+1}: Expected m/z (M+ion)+: {parent_mass:.4f}"
                )
            elif self.mode == "Negative":
                parent_mass = self.compound_MW - MW
                self.logger.info(
                    f"Ion {ion_index+1}: Expected m/z (M-ion)-: {parent_mass:.4f}"
                )

            # Generate test masses based on tolerance
            test_ion_masses = [
                parent_mass - tolerance,
                parent_mass,
                parent_mass + tolerance,
            ]
            self.logger.info(
                f"Testing m/z values: {[round(m) for m in test_ion_masses]}"
            )

            # Try to find exact matches first
            ion_matches = []
            for ion_mass in test_ion_masses:
                ion_match_max = self.get_max_match_indices(ion_mass)
                ion_matches.append(ion_match_max)

            # If no exact matches found, try searching in any position
            if not any([ion_match[1] for ion_match in ion_matches]):
                self.logger.info("No maximum matches found, searching for any matches")
                for ion_mass in test_ion_masses:
                    ion_matches = []
                    ion_any_match = self.get_any_match_indices(ion_mass)
                    ion_matches.append(ion_any_match)

            # Process matches if found
            if any([ion_max_match[1] for ion_max_match in ion_matches]):
                self.logger.info(
                    f"Found matches for compound with m/z ~{round(parent_mass)}"
                )
                for ion_match in ion_matches:
                    if ion_match[1]:
                        self.logger.info(
                            f"Match at m/z {ion_match[0]}, indices: {ion_match[1]}"
                        )
                        EIC_data.append(
                            self.get_extracted_ion_count(
                                RT_data=self.MSdata["RT"],
                                mz_data=self.MSdata["mz_data"],
                                ion_mass=ion_match[0],
                            )
                        )
                        max_mz_match.append(ion_match[2])
                        ions_matched.append(
                            (get_smiles(ion_to_alter_mol), ion_match[0])
                        )
                        RT_matched.append(self.MSpeakdata["RT"][ion_match[1]])
                        TIC_matched.append(self.MSpeakdata["TIC"][ion_match[1]])
                        mz_data_matched.append(
                            [
                                self.MSpeakdata["mz_data"][index]
                                for index in ion_match[1]
                            ]
                        )
                        mz_strongest_value.append(
                            self.get_strongest_mz_pattern(
                                [
                                    self.MSpeakdata["mz_data"][index]
                                    for index in ion_match[1]
                                ]
                            )
                        )
            else:
                self.logger.info(f"No matches found for m/z {round(parent_mass)}")

        # Calculate maximum EIC signal if data exists
        if EIC_data:
            max_EIC_signal = np.max([data[1] for EIC in EIC_data for data in EIC])
            self.logger.info(f"Maximum EIC signal: {max_EIC_signal:.2f}")
        else:
            max_EIC_signal = None
            self.logger.warning("No EIC data found for this compound")

        # Store analysis results
        self.analysedata = {
            "EIC_data": EIC_data,
            "max_EIC_signal": max_EIC_signal,
            "max_mz_match": max_mz_match,
            "ions": ions_matched,
            "RT": RT_matched,
            "TIC": TIC_matched,
            "mz_data": mz_data_matched,
            "mz_strongest": mz_strongest_value,
        }

        matched_ions_str = (
            ", ".join([f"{ion[0]}({ion[1]})" for ion in ions_matched])
            if ions_matched
            else "None"
        )
        self.logger.info(f"Analysis complete. Matched ions: {matched_ions_str}")

        # Find RT values for maximum TIC signals
        if RT_matched and TIC_matched:
            self.max_RT_ion_match()

        return self.analysedata

    def get_strongest_mz_pattern(self, mz_values: list) -> largest_mz_pattern:
        """
        Finds largest mz signal from list of mz patterns and returns mz max value,
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

    def calculate_eic_area(self) -> float:
        """Calculate the area under the most dominant extracted ion count curve using the trapezoidal rule"""
        self.logger.info("Calculating EIC area")

        eic_data = self.analysedata.get("EIC_data", [])
        if not eic_data:
            self.logger.warning("No EIC data available for area calculation")
            return 0.0

        # Find the most dominant EIC curve (the one with the highest peak)
        max_peak_height = 0
        dominant_eic_index = 0

        for i, eic in enumerate(eic_data):
            if not eic:
                continue

            peak_height = max(point[1] for point in eic)
            if peak_height > max_peak_height:
                max_peak_height = peak_height
                dominant_eic_index = i

        # If we found a dominant curve
        if max_peak_height > 0 and dominant_eic_index < len(eic_data):
            dominant_curve = eic_data[dominant_eic_index]
            self.logger.info(
                f"Using dominant EIC curve at index {dominant_eic_index} with peak height {max_peak_height:.2f}"
            )

            # Sort by retention time to ensure correct integration
            dominant_curve_sorted = sorted(dominant_curve, key=lambda point: point[0])

            # Calculate area using trapezoidal rule
            area = 0.0
            for i in range(1, len(dominant_curve_sorted)):
                rt_diff = dominant_curve_sorted[i][0] - dominant_curve_sorted[i - 1][0]
                height_avg = (
                    dominant_curve_sorted[i][1] + dominant_curve_sorted[i - 1][1]
                ) / 2
                area += rt_diff * height_avg

            self.logger.info(f"Calculated EIC area: {area:.2f}")
            # Update analysedata with the calculated area
            self.analysedata["EIC_area"] = area

            return area
        else:
            self.logger.warning("No dominant EIC curve found")
            self.analysedata["EIC_area"] = 0.0
            return 0.0

    def max_RT_ion_match(self) -> dict:
        """
        Find retention times corresponding to maximum TIC signals for each matched ion.
        
        Returns:
            dict: Dictionary containing RT values of maximum intensity for each matched ion,
                with the overall maximum RT and its corresponding ion
        """
        self.logger.info("Finding RT values for maximum TIC signals")
        
        # Check if we have RT and TIC data from previous analysis
        if not self.analysedata or "RT" not in self.analysedata or not self.analysedata["RT"]:
            self.logger.warning("No RT data available from analysis")
            return {
                "rt_values": [],
                "max_rt": None,
                "max_ion": None
            }
        
        rt_at_max = []
        max_intensities = []
        ions = []
        
        # Process each RT array and corresponding TIC array
        for i, (rt_array, tic_array) in enumerate(zip(self.analysedata["RT"], self.analysedata["TIC"])):
            if len(rt_array) > 0 and len(tic_array) > 0:
                # Find index of maximum TIC value
                max_idx = np.argmax(tic_array)
                # Get corresponding RT value
                rt_max = rt_array[max_idx]
                # Get intensity at maximum
                intensity_max = tic_array[max_idx]
                # Get ion information if available
                ion = self.analysedata["ions"][i] if i < len(self.analysedata["ions"]) else None
                
                rt_at_max.append(rt_max)
                max_intensities.append(intensity_max)
                ions.append(ion)
                
                self.logger.debug(f"Ion match {i+1}: RT at max TIC = {rt_max:.2f} min, Max TIC = {intensity_max:.2f}")
        
        # Find the overall maximum
        overall_max_rt = None
        overall_max_ion = None
        
        if max_intensities:
            overall_max_idx = np.argmax(max_intensities)
            overall_max_rt = rt_at_max[overall_max_idx]
            overall_max_ion = ions[overall_max_idx]
            
            self.logger.info(f"Overall maximum TIC at RT {overall_max_rt:.2f} min for ion {overall_max_ion}")
        
        # Store results in analysedata for easy access
        self.analysedata["rt_at_max_tic"] = rt_at_max
        self.analysedata["max_rt"] = overall_max_rt
        self.analysedata["max_ion"] = overall_max_ion
        
        return {
            "rt_values": rt_at_max,
            "max_rt": overall_max_rt,
            "max_ion": overall_max_ion
        }

    def create_report(self, folder: str = "reports", compound_name: str = None) -> str:
        """Create a report for the analyzed spectrum"""
        if not compound_name:
            compound_name = get_path_leaf(self._filepath)

        # Create a report generator instance
        report_generator = MSReport(output_dir=folder)

        # Generate the report
        return report_generator.create_compound_report(
            msmode=self.mode,
            RT_values=self.MSdata["RT"],
            TIC_values=self.MSdata["TIC"],
            compound_name=compound_name,
            mol=self.compound_mol,
            analysedata=self.analysedata,
        )
