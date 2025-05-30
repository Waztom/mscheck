from datetime import datetime
import os
import re
import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
import yaml
from rdkit.Chem import Descriptors, AllChem


# Local imports
from analyse import AnalyseSpectrum
from report import MSReport
from heatmap import MSHeatmapGenerator
from logging_config import setup_logger, get_logger


class BulkAnalyser:
    """Handles bulk analysis of mass spectrometry data using CSV input"""

    def __init__(self, config_path: str):
        """
        Initialize bulk analyzer with a configuration file

        Args:
            config_path: Path to YAML configuration file
        """
        # Load the configuration
        self.config = self.load_config(config_path)

        # Extract paths from config
        self.csv_input_path = self.config["paths"]["csv_file"]
        self.data_dir = self.config["paths"].get(
            "data_dir", os.path.dirname(self.csv_input_path)
        )
        self.report_dir = (
            self.config["paths"]
            .get("output", {})
            .get(
                "base_dir",
                os.path.join(os.path.dirname(self.csv_input_path), "reports"),
            )
        )

        self.batch_data = None
        self.processed_samples = 0
        self.total_samples = 0
        self.errors = []
        self.logger = logging.getLogger("BulkAnalyser")

        self.logger.info(f"Initialized BulkAnalyser with config: {config_path}")
        self.logger.info(f"CSV file: {self.csv_input_path}")
        self.logger.info(f"Data directory: {self.data_dir}")
        self.logger.info(f"Report directory: {self.report_dir}")

    def load_config(self, config_path: str) -> Dict:
        """
        Load configuration file

        Args:
            config_path: Path to YAML configuration file

        Returns:
            Dictionary with configuration

        Raises:
            FileNotFoundError: If config file doesn't exist
            ValueError: If config is invalid or missing required fields
        """
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        try:
            with open(config_path, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)

            # Validate required sections
            if "paths" not in config:
                raise ValueError("Configuration missing 'paths' section")

            if "csv_file" not in config["paths"]:
                raise ValueError("Configuration missing 'paths.csv_file' setting")

            # Set default parameters if not provided
            if "parameters" not in config:
                config["parameters"] = {}

            if "analysis_types" not in config["parameters"]:
                config["parameters"]["analysis_types"] = [
                    "product",
                    "reactant",
                    "internal-std",
                ]

            if "modes" not in config["parameters"]:
                config["parameters"]["modes"] = ["Positive"]

            if "tolerance" not in config["parameters"]:
                config["parameters"]["tolerance"] = 1

            return config

        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing YAML configuration: {str(e)}")
        except Exception as e:
            raise ValueError(f"Error loading configuration: {str(e)}")

    def load_data(self) -> pd.DataFrame:
        """
        Load data from CSV file

        Returns:
            DataFrame with loaded data
        """
        try:
            self.batch_data = pd.read_csv(self.csv_input_path)
            self.total_samples = len(self.batch_data)
            self.logger.info(
                f"Successfully loaded {self.total_samples} samples from {self.csv_input_path}"
            )
            return self.batch_data
        except Exception as e:
            error_msg = f"Error loading CSV file: {str(e)}"
            self.errors.append(error_msg)
            self.logger.error(error_msg)
            self.batch_data = pd.DataFrame()
            return self.batch_data

    def analyse_batch(
        self,
        analysis_types: List[str] = ["product", "intermediate", "reactant"],
        modes: List[str] = ["Positive"],
        tolerance: int = 1,
    ) -> pd.DataFrame:
        """
        Process all samples in the batch

        Args:
            analysis_types: List of analysis types to perform (e.g., product, reactant)
            modes: List of ionization modes to analyze
            tolerance: Default mass tolerance if not specified in CSV

        Returns:
            DataFrame with analysis results
        """
        if self.batch_data is None:
            self.load_data()

        if len(self.batch_data) == 0:
            self.logger.warning("No data to analyze")
            return self.batch_data

        # Create report directories
        for mode in modes:
            os.makedirs(os.path.join(self.report_dir, mode), exist_ok=True)

        # Track processing statistics
        self.processed_samples = 0
        successful_analyses = 0
        total_analyses_attempted = 0

        # Process each row in the CSV
        for batch_index, batch_row in self.batch_data.iterrows():
            try:
                sample_id = str(batch_row.get("sample-ID", f"sample_{batch_index}"))
                mzML_filename = batch_row["mzML-filename"]
                self.logger.info(
                    f"Processing sample {sample_id} ({batch_index + 1}/{len(self.batch_data)})"
                )

                # Find mzML file
                mzML_filepath = self.find_mzml_file(mzML_filename)
                if mzML_filepath is None:
                    error_msg = f"Warning: mzML file not found for sample {sample_id}: {mzML_filename}"
                    self.errors.append(error_msg)
                    self.logger.warning(error_msg)
                    continue

                # Process each analysis type
                for analysis_type in analysis_types:
                    try:
                        # Get number of compounds of this type by scanning columns
                        no_compounds = self.get_compound_count(batch_row, analysis_type)
                        if no_compounds <= 0:
                            continue

                        # Get ion information
                        ions_add_column = f"{analysis_type}-ions-to-add"
                        ions_sub_column = f"{analysis_type}-ions-to-sub"

                        # Get ions with safe parsing
                        analysis_type_ions_to_add = self.parse_ions(
                            batch_row.get(ions_add_column, "")
                        )
                        analysis_type_ions_to_sub = self.parse_ions(
                            batch_row.get(ions_sub_column, "")
                        )

                        # Get tolerance - use default if not specified
                        tolerance_column = f"{analysis_type}-match-tolerance"
                        if tolerance_column in batch_row:
                            analysis_type_match_tolerance = int(
                                batch_row[tolerance_column]
                            )
                        else:
                            analysis_type_match_tolerance = tolerance

                        # Process each compound of this type
                        for compound_idx in range(no_compounds):
                            # Get compound SMILES
                            compound_column = f"{analysis_type}-{compound_idx + 1}"
                            if compound_column not in batch_row:
                                self.logger.warning(
                                    f"'{compound_column}' column missing for sample {sample_id}"
                                )
                                continue

                            analysis_type_smiles = batch_row[compound_column]

                            # Add this after getting analysis_type_smiles
                            if analysis_type == "internal-std":
                                # For internal standards, look up SMILES from config
                                std_name = analysis_type_smiles  # The value is actually the name
                                for is_info in (
                                    self.config.get("conversion", {})
                                    .get("internal_standards", {})
                                    .get("standards", [])
                                ):
                                    if is_info.get("name") == std_name:
                                        analysis_type_smiles = is_info.get("smiles")
                                        self.logger.info(
                                            f"Found SMILES for {std_name}: {analysis_type_smiles}"
                                        )
                                        break

                            if (
                                not isinstance(analysis_type_smiles, str)
                                or not analysis_type_smiles.strip()
                            ):
                                self.logger.warning(
                                    f"Empty SMILES in {compound_column} for sample {sample_id}"
                                )
                                continue

                            # Desalt SMILES before analysis
                            original_smiles = analysis_type_smiles
                            desalted_smiles = self.desalt_smiles(original_smiles)
                            if desalted_smiles != original_smiles:
                                self.logger.info(
                                    f"Sample {sample_id}: Desalted {compound_column} from {original_smiles} → {desalted_smiles}"
                                )
                                analysis_type_smiles = desalted_smiles

                            # Log the original SMILES
                            self.logger.info(
                                f"Analyzing {analysis_type}-{compound_idx+1}: {analysis_type_smiles}"
                            )

                            # Calculate and log molecular weight
                            try:
                                mol = Chem.MolFromSmiles(analysis_type_smiles)
                                if mol:
                                    mw = Descriptors.MolWt(mol)
                                    formula = AllChem.CalcMolFormula(mol)
                                    self.logger.info(
                                        f"Molecular weight: {mw:.2f} Da, Formula: {formula}"
                                    )

                                    # Calculate m/z values to search for (based on mode)
                                    if mode == "Positive":
                                        mz_protonated = mw + 1.007825  # Mass of proton
                                        self.logger.info(
                                            f"Expected m/z (M+H)+: {mz_protonated:.4f}"
                                        )
                                    else:  # Negative mode
                                        mz_deprotonated = (
                                            mw - 1.007825
                                        )  # Loss of proton
                                        self.logger.info(
                                            f"Expected m/z (M-H)-: {mz_deprotonated:.4f}"
                                        )
                            except Exception as e:
                                self.logger.warning(
                                    f"Error calculating MW for {analysis_type_smiles}: {str(e)}"
                                )

                            # Process each ionization mode
                            for mode in modes:
                                try:
                                    total_analyses_attempted += 1
                                    self.logger.info(
                                        f"Analyzing {analysis_type}-{compound_idx + 1} in {mode} mode"
                                    )

                                    # Perform analysis
                                    analysis_obj = AnalyseSpectrum(
                                        mzMLfilepath=mzML_filepath, mode=mode
                                    )

                                    # Set and log default ions if none provided, based on mode
                                    if not analysis_type_ions_to_add:
                                        analysis_type_ions_to_add = ["[H]"]
                                        self.logger.info(
                                            f"No ions specified, using default for {mode} mode: {analysis_type_ions_to_add}"
                                        )
                                    else:
                                        self.logger.info(
                                            f"Using specified ions for {mode} mode: {analysis_type_ions_to_add}"
                                        )

                                    # Log actual SMILES representation of ions being addeed
                                    self.logger.info(
                                        f"Final ion SMILES for {mode} mode: {analysis_type_ions_to_add}"
                                    )

                                    # Now perform the analysis
                                    analysis_obj.analyse(
                                        compoundsmiles=analysis_type_smiles,
                                        ionstoadd=analysis_type_ions_to_add,
                                        ionstosub=analysis_type_ions_to_sub,
                                        tolerance=analysis_type_match_tolerance,
                                    )

                                    # Calculate EIC area
                                    eic_area = analysis_obj.calculate_eic_area()
                                    self.logger.info(f"EIC Area: {eic_area:.2f}")

                                    # Create report
                                    report_subdir = os.path.join(self.report_dir, mode)
                                    report_name = f"{sample_id}_{analysis_type}_{compound_idx + 1}"
                                    analysis_obj.create_report(
                                        folder=report_subdir, compound_name=report_name
                                    )

                                    # Update results in DataFrame using proper column names
                                    signal_column = f"{analysis_type}-{compound_idx + 1}-max-EIC-signal-{mode}"
                                    mz_match_column = f"{analysis_type}-{compound_idx + 1}-max-mz-match-{mode}"
                                    ions_column = f"{analysis_type}-{compound_idx + 1}-ions-matched-{mode}"
                                    eic_area_column = f"{analysis_type}-{compound_idx + 1}-EIC-area-{mode}"

                                    # Use loc to properly update DataFrame
                                    self.batch_data.loc[batch_index, signal_column] = (
                                        analysis_obj.analysedata["max_EIC_signal"]
                                    )
                                    self.batch_data.loc[
                                        batch_index, mz_match_column
                                    ] = str(analysis_obj.analysedata["max_mz_match"])
                                    self.batch_data.loc[batch_index, ions_column] = str(
                                        analysis_obj.analysedata["ions"]
                                    )
                                    self.batch_data.loc[
                                        batch_index, eic_area_column
                                    ] = eic_area

                                    successful_analyses += 1
                                except Exception as e:
                                    error_msg = f"Error analyzing {analysis_type}-{compound_idx + 1} in {mode} for sample {sample_id}: {str(e)}"
                                    self.errors.append(error_msg)
                                    self.logger.error(error_msg)
                    except Exception as e:
                        error_msg = f"Error processing {analysis_type} for sample {sample_id}: {str(e)}"
                        self.errors.append(error_msg)
                        self.logger.error(error_msg)

                self.processed_samples += 1

            except Exception as e:
                error_msg = f"Error processing sample {batch_index}: {str(e)}"
                self.errors.append(error_msg)
                self.logger.error(error_msg)

        # Print summary
        self.logger.info(
            f"Analysis complete: {successful_analyses}/{total_analyses_attempted} analyses successful"
        )
        self.logger.info(
            f"Processed {self.processed_samples}/{self.total_samples} samples"
        )

        if self.errors:
            self.logger.warning(
                f"Encountered {len(self.errors)} errors during processing"
            )

        return self.batch_data

    def find_mzml_file(
        self, filename: str, custom_data_dir: str = None
    ) -> Optional[str]:
        """
        Helper method to locate mzML files with flexible path handling

        Args:
            filename: Base filename with or without extension
            custom_data_dir: Optional custom directory to search in (defaults to self.data_dir)

        Returns:
            Full path to mzML file if found, None otherwise
        """
        # Use custom directory if provided, otherwise use default
        data_dir = custom_data_dir if custom_data_dir is not None else self.data_dir

        # Try different possible file paths
        possible_paths = [
            os.path.join(data_dir, f"{filename}.mzML"),
            os.path.join(data_dir, filename),
        ]

        self.logger.info(f"Searching for mzML file: {filename}")
        self.logger.info(f"Looking in directory: {data_dir}")
        self.logger.info(f"Possible paths: {possible_paths}")

        # Try with extension if not already present
        if not filename.endswith(".mzML"):
            possible_paths.append(os.path.join(data_dir, f"{filename}.mzML"))
        else:
            # If filename already has extension, try without it
            base_name = filename[:-5]
            possible_paths.append(os.path.join(data_dir, base_name))

        # Check if any path exists
        for path in possible_paths:
            if os.path.exists(path):
                self.logger.info(f"Found mzML file at: {path}")
                return path

        self.logger.warning(f"Could not find mzML file for {filename}")
        return None

    def parse_ions(self, ion_str: str) -> List[str]:
        """
        Helper method to parse ion strings from CSV

        Args:
            ion_str: String containing ions (may be H, empty, or comma-separated)

        Returns:
            List of ion SMILES strings
        """
        if not isinstance(ion_str, str):
            return []

        ion_str = ion_str.strip()

        # Handle special case for hydrogen
        if ion_str == "H":
            return ["[H]"]
        elif not ion_str:
            return []

        # Handle the case where the string starts with the ion without brackets
        if ion_str.startswith("[") and "]" in ion_str:
            # Already formatted ions
            return [s.strip() for s in ion_str.split(",") if s.strip()]
        else:
            # Need to split and add brackets if needed
            ions = []
            for ion in ion_str.split(","):
                ion = ion.strip()
                if ion:
                    # Add brackets if needed
                    if not (ion.startswith("[") and ion.endswith("]")):
                        ion = f"[{ion}]"
                    ions.append(ion)
            return ions

    def save_results(self, output_path: str) -> bool:
        """
        Save analysis results to CSV

        Args:
            output_path: Path to save CSV results

        Returns:
            Boolean indicating success
        """
        if self.batch_data is None or len(self.batch_data) == 0:
            self.logger.warning("No data to save")
            return False

        try:
            # Handle directory vs file path
            if os.path.isdir(output_path) or not output_path.lower().endswith(".csv"):
                os.makedirs(output_path, exist_ok=True)
                timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
                output_file = os.path.join(
                    output_path, f"analysis_results_{timestamp}.csv"
                )
            else:
                output_file = output_path
                os.makedirs(os.path.dirname(output_file), exist_ok=True)

            self.batch_data.to_csv(output_file, index=False)
            self.logger.info(f"Results saved to {output_file}")
            return True

        except Exception as e:
            error_msg = f"Error saving results: {str(e)}"
            self.errors.append(error_msg)
            self.logger.error(error_msg)
            return False

    def get_summary(self) -> Dict:
        """
        Get summary statistics from the analysis

        Returns:
            Dictionary with summary information
        """
        return {
            "total_samples": self.total_samples,
            "processed_samples": self.processed_samples,
            "errors": len(self.errors),
            "error_messages": self.errors,
        }

    def visualize_results(
        self, metrics: list = None, output_dir: str = None, show_plots: bool = False
    ) -> None:
        """Generate visualizations for analysis results"""
        if output_dir is None:
            output_dir = self.report_dir

        # Create visualization generator
        vis = MSHeatmapGenerator(data=self.batch_data, output_dir=output_dir)

        # Generate all visualizations
        vis.generate_all_visualizations(metrics=metrics, show_plots=show_plots)

        # Add any errors to the bulk analyzer's error list
        self.errors.extend(vis.errors)

    def process_samples(
        self,
        analysis_types: List[str] = ["product", "reactant", "internal-std"],
        modes: List[str] = ["Positive"],
        tolerance: int = 1,
    ):
        """
        Process all samples and generate individual reports

        Args:
            analysis_types: List of compound types to analyze (e.g., product, reactant)
            modes: List of ionization modes (e.g., Positive, Negative)
            tolerance: Default mass tolerance in ppm

        Returns:
            DataFrame with analysis results
        """
        if self.batch_data is None:
            self.load_data()

        # Process the samples with the specified parameters
        self.analyse_batch(
            analysis_types=analysis_types, modes=modes, tolerance=tolerance
        )

        # Save the results
        results_path = os.path.join(
            os.path.dirname(self.csv_input_path), "analysis_results.csv"
        )
        self.save_results(results_path)

        # Print analysis summary
        summary = self.get_summary()
        self.logger.info("\nAnalysis Summary:")
        self.logger.info(f"Total samples: {summary['total_samples']}")
        self.logger.info(f"Processed samples: {summary['processed_samples']}")
        self.logger.info(f"Errors encountered: {summary['errors']}")

        return self.batch_data

    def generate_visualizations(self):
        """Generate heatmaps and other visualizations"""
        self.logger.info("\nGenerating heatmaps and visualizations...")
        heatmap_generator = MSHeatmapGenerator(
            data=self.batch_data, output_dir=self.report_dir
        )
        heatmap_generator.generate_all_visualizations(show_plots=False)
        self.errors.extend(heatmap_generator.errors)

    def extract_sample_data(self):
        """Extract organized sample data from batch results"""
        self.logger.info("Extracting sample data for reports...")
        samples = {}

        for batch_index, row in self.batch_data.iterrows():
            sample_id = row.get("sample-ID", "unknown")
            if sample_id not in samples:
                samples[sample_id] = {
                    "compounds": [],
                    "RT_values": None,
                    "TIC_values": None,
                    "all_mz_data": None,
                }

            # Get sample chromatogram data if not already loaded
            self._load_sample_chromatogram(samples, sample_id, row)

            # Process each compound type
            for compound_type in [
                "reactant",
                "product",
                "internal-std",
                "intermediate",
            ]:
                self._extract_compound_data(samples, sample_id, row, compound_type)

        self.logger.info(f"Found {len(samples)} unique samples with compounds")
        return samples

    def _load_sample_chromatogram(self, samples, sample_id, row):
        """Helper method to load chromatogram data for a sample"""
        if samples[sample_id]["RT_values"] is not None:
            return  # Already loaded

        try:
            # Find mzML file
            mzML_filename = row["mzML-filename"]
            mzML_filepath = self.find_mzml_file(mzML_filename)
            if mzML_filepath is not None:
                # Create analyzer to get RT/TIC
                mode = "Positive"
                temp_analyzer = AnalyseSpectrum(mzMLfilepath=mzML_filepath, mode=mode)
                samples[sample_id]["RT_values"] = temp_analyzer.MSdata["RT"]
                samples[sample_id]["TIC_values"] = temp_analyzer.MSdata["TIC"]

                # Store MZ data if available
                if "mz_data" in temp_analyzer.MSdata:
                    samples[sample_id]["all_mz_data"] = temp_analyzer.MSdata["mz_data"]
                    self.logger.info(f"Added complete MZ data to sample {sample_id}")
        except Exception as e:
            self.logger.error(f"Error loading RT/TIC for sample {sample_id}: {str(e)}")

    def _extract_compound_data(self, samples, sample_id, row, compound_type):
        """Extract compound data for a specific compound type"""
        # Get number of compounds by scanning column headers
        no_compounds = self.get_compound_count(row, compound_type)

        # Process each compound
        for compound_idx in range(no_compounds):
            compound_data = self._extract_single_compound(
                sample_id, row, compound_type, compound_idx
            )

            if compound_data:
                samples[sample_id]["compounds"].append(compound_data)

    def _extract_single_compound(self, sample_id, row, compound_type, compound_idx):
        """Extract data for a single compound"""
        # Get compound information
        compound_column = f"{compound_type}-{compound_idx + 1}"
        smiles_value = row.get(compound_column)

        # Skip if no SMILES
        if not isinstance(smiles_value, str) or not smiles_value.strip():
            return None

        # Desalt SMILES
        smiles_value = self.desalt_smiles(smiles_value)

        # Get analysis results for Positive mode
        mode = "Positive"
        signal_column = f"{compound_type}-{compound_idx + 1}-max-EIC-signal-{mode}"
        eic_area_column = f"{compound_type}-{compound_idx + 1}-EIC-area-{mode}"
        mz_match_column = f"{compound_type}-{compound_idx + 1}-max-mz-match-{mode}"

        # Skip if analysis wasn't completed
        if signal_column not in row or pd.isna(row[signal_column]):
            return None

        # Create RDKit molecule
        try:
            mol = Chem.MolFromSmiles(smiles_value)
            if mol is None:
                self.logger.warning(
                    f"Could not create molecule from SMILES: {smiles_value}"
                )
                return None
        except Exception as e:
            self.logger.error(f"Error creating molecule: {str(e)}")
            return None

        # Get the compound's mass spectrum and other data
        rt_max, intensity_max, eic_data, mz_strongest, all_mz_data = (
            self._get_compound_ms_data(
                sample_id, row, compound_type, compound_idx, smiles_value, mode
            )
        )

        # Create compound entry
        return {
            "compound_name": f"{compound_type.capitalize()} {compound_idx + 1}",
            "compound_type": compound_type,
            "msmode": mode,
            "mol": mol,
            "smiles": smiles_value,
            "rt_max": rt_max,
            "intensity_max": intensity_max,
            "EIC_data": eic_data,
            "signal": row.get(signal_column),
            "EIC_area": row.get(eic_area_column),
            "max_mz_match": row.get(mz_match_column, "False") == "True",
            "mz_strongest": mz_strongest,
            "all_mz_data": all_mz_data,
        }

    def _get_compound_ms_data(
        self, sample_id, row, compound_type, compound_idx, smiles_value, mode
    ):
        """Get mass spectrum data for a compound"""
        rt_max = None
        intensity_max = None
        eic_data = None
        mz_strongest = None
        all_mz_data = None

        # Construct report path
        report_subdir = os.path.join(self.report_dir, mode)
        report_name = f"{sample_id}_{compound_type}_{compound_idx + 1}"
        report_path = os.path.join(report_subdir, f"{report_name}-report.svg")

        # Only proceed if report exists
        if not os.path.exists(report_path):
            return rt_max, intensity_max, eic_data, mz_strongest, all_mz_data

        try:
            # Get mzML file
            mzML_filename = row["mzML-filename"]
            mzML_filepath = self.find_mzml_file(mzML_filename)

            # Create analyzer and run analysis
            temp_analyzer = AnalyseSpectrum(mzMLfilepath=mzML_filepath, mode=mode)

            # Get parameters for analysis
            ionstoadd = self.parse_ions(row.get(f"{compound_type}-ions-to-add", ""))
            ionstosub = self.parse_ions(row.get(f"{compound_type}-ions-to-sub", ""))
            tolerance = int(row.get(f"{compound_type}-match-tolerance", 1))

            # Run analysis
            temp_analyzer.analyse(
                compoundsmiles=smiles_value,
                ionstoadd=ionstoadd,
                ionstosub=ionstosub,
                tolerance=tolerance,
            )

            # Extract RT max and intensity max
            if (
                "RT" in temp_analyzer.analysedata
                and len(temp_analyzer.analysedata["RT"]) > 0
            ):
                rt_values = temp_analyzer.analysedata["RT"][0]
                tic_values = temp_analyzer.analysedata["TIC"][0]
                if len(rt_values) > 0 and len(tic_values) > 0:
                    max_idx = np.argmax(tic_values)
                    rt_max = rt_values[max_idx]
                    intensity_max = tic_values[max_idx]

            # Extract EIC data
            if (
                "EIC_data" in temp_analyzer.analysedata
                and len(temp_analyzer.analysedata["EIC_data"]) > 0
            ):
                eic_data = temp_analyzer.analysedata["EIC_data"][0]

                # Add this to log details about matches
                if "matched_ions" in temp_analyzer.analysedata:
                    matches = temp_analyzer.analysedata["matched_ions"]
                    for i, match in enumerate(matches):
                        formula = match.get("formula", "Unknown")
                        mz_theoretical = match.get("mz_theoretical", 0)
                        mz_found = match.get("mz_found", 0)
                        ppm_error = match.get("ppm_error", 0)
                        adduct = match.get("adduct", "Unknown")

                        self.logger.info(f"Match {i+1}: {formula} as {adduct}")
                        self.logger.info(
                            f"  m/z theoretical: {mz_theoretical:.4f}, found: {mz_found:.4f}, error: {ppm_error:.2f} ppm"
                        )

            # Extract mass spectrum data
            if (
                "mz_strongest" in temp_analyzer.analysedata
                and len(temp_analyzer.analysedata["mz_strongest"]) > 0
            ):
                mz_masses, mz_intensities, _ = temp_analyzer.analysedata[
                    "mz_strongest"
                ][0]
                mz_strongest = (mz_masses, mz_intensities)

            # Get full MZ data
            if hasattr(temp_analyzer, "MSdata") and "mz_data" in temp_analyzer.MSdata:
                all_mz_data = temp_analyzer.MSdata["mz_data"]

        except Exception as e:
            self.logger.error(f"Error extracting MS data: {str(e)}")

        return rt_max, intensity_max, eic_data, mz_strongest, all_mz_data

    def generate_compound_reports(self, samples=None):
        """Generate multi-compound reports for samples"""
        if samples is None:
            samples = self.extract_sample_data()

        report_generator = MSReport(
            output_dir=os.path.join(self.report_dir, "compound_reports")
        )
        report_paths = {}

        # First pass - collect information about all reports we'll generate
        available_reports = []
        for sample_id, sample_data in samples.items():
            if (
                len(sample_data.get("compounds", [])) > 0
                and sample_data.get("RT_values") is not None
            ):
                report_filename = f"Sample_{sample_id}_Analysis.html"
                report_path = os.path.join(
                    self.report_dir, "compound_reports", f"{sample_id}_analysis.html"
                )
                report_title = f"Sample {sample_id} Analysis"

                available_reports.append(
                    {
                        "title": report_title,
                        "path": report_filename,
                        "sample_id": sample_id,
                    }
                )

        # Second pass - actually generate the reports with navigation
        for i, report_info in enumerate(available_reports):
            sample_id = report_info["sample_id"]
            sample_data = samples[sample_id]

            # Filter valid compounds
            valid_compounds = [
                c
                for c in sample_data["compounds"]
                if c.get("rt_max") is not None and c.get("mol") is not None
            ]

            if not valid_compounds:
                self.logger.warning(
                    f"No compounds with valid RT data for sample {sample_id}"
                )
                continue

            self.logger.info(
                f"Creating report with {len(valid_compounds)} compounds for {sample_id}"
            )

            # Generate report
            report_path = report_generator.create_annotated_tic_report(
                RT_values=sample_data["RT_values"],
                TIC_values=sample_data["TIC_values"],
                compounds=valid_compounds,
                report_title=report_info["title"],
                html_output=True,
                svg_layout="triple",
                available_reports=available_reports,
                current_report_index=i,
            )

            if report_path:
                report_paths[sample_id] = report_path
                self.logger.info(f"Generated report: {report_path}")
            else:
                self.logger.error(f"Failed to generate report for sample {sample_id}")

        return report_paths

    def generate_conversion_reports(self, conversion_df):
        """
        Generate reports and visualizations for conversion results

        Args:
            conversion_df: DataFrame with conversion results

        Returns:
            Dictionary of report paths
        """
        self.logger.info("Generating conversion reports...")

        # Check if we have results
        if conversion_df is None or len(conversion_df) == 0:
            self.logger.warning("No conversion data available for reports")
            return {}

        # Get report settings
        report_dir = os.path.join(self.report_dir, "conversion_reports")
        self.logger.info(f"Creating conversion reports directory: {report_dir}")

        # Create report directory
        os.makedirs(report_dir, exist_ok=True)

        # Create unique timestamp for this report set
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_paths = {}

        # ===== 1. Generate detailed CSV summary =====
        summary_path = os.path.join(report_dir, f"conversion_summary_{timestamp}.csv")
        conversion_df.to_csv(summary_path, index=False)
        report_paths["detailed_summary"] = summary_path
        self.logger.info(f"Saved detailed conversion summary to: {summary_path}")

        # ===== 2. Generate summary by sample =====
        try:
            # Group by sample and calculate statistics
            sample_summary = (
                conversion_df.groupby("sample_id")
                .agg(
                    {
                        "concentration_reactant_uM": [
                            "mean",
                            "std",
                            "min",
                            "max",
                            "count",
                        ]
                    }
                )
                .reset_index()
            )

            # Flatten column names
            sample_summary.columns = [
                "sample_id",
                "mean_conversion",
                "std_conversion",
                "min_conversion",
                "max_conversion",
                "measurement_count",
            ]

            sample_path = os.path.join(report_dir, f"sample_summary_{timestamp}.csv")
            sample_summary.to_csv(sample_path, index=False)
            report_paths["sample_summary"] = sample_path
            self.logger.info(f"Saved sample conversion summary to: {sample_path}")
        except Exception as e:
            self.logger.error(f"Error generating sample summary: {str(e)}")

        # ===== 3. Generate summary by reactant =====
        try:
            # Group by reactant and calculate statistics
            reactant_summary = (
                conversion_df.groupby("reactant_name")
                .agg({"conversion_percent": ["mean", "std", "min", "max", "count"]})
                .reset_index()
            )

            # Flatten column names
            reactant_summary.columns = [
                "reactant_name",
                "mean_conversion",
                "std_conversion",
                "min_conversion",
                "max_conversion",
                "measurement_count",
            ]

            reactant_path = os.path.join(
                report_dir, f"reactant_summary_{timestamp}.csv"
            )
            reactant_summary.to_csv(reactant_path, index=False)
            report_paths["reactant_summary"] = reactant_path
            self.logger.info(f"Saved reactant conversion summary to: {reactant_path}")
        except Exception as e:
            self.logger.error(f"Error generating reactant summary: {str(e)}")

        self.logger.info(
            f"Generated {len(report_paths)} conversion reports in {report_dir}"
        )
        return report_paths

    def run_complete_workflow(self):
        """
        Run the complete analysis workflow using parameters from config

        Returns:
            Dictionary mapping sample IDs to their report paths
        """
        # Get parameters from config
        analysis_types = self.config["parameters"]["analysis_types"]
        modes = self.config["parameters"]["modes"]
        tolerance = self.config["parameters"]["tolerance"]

        self.logger.info("Starting complete analysis workflow...")
        self.logger.info(f"Analysis types: {', '.join(analysis_types)}")
        self.logger.info(f"Modes: {', '.join(modes)}")
        self.logger.info(f"Default tolerance: {tolerance}")

        # 1. Process all samples with the specified parameters
        self.process_samples(
            analysis_types=analysis_types, modes=modes, tolerance=tolerance
        )

        # 2. Calculate response factors and conversions if enabled
        if "conversion" in self.config and self.config["conversion"].get(
            "enabled", True
        ):
            self.logger.info("Conversion analysis is enabled in config")
            self.logger.info("Calculating response factors for conversion analysis...")

            response_factors = self.calculate_response_factors()
            self.logger.info(f"Found {len(response_factors)} response factors")

            if response_factors:
                self.logger.info("Calculating conversions using response factors...")
                conversion_results = self.calculate_conversions(response_factors)

                if not conversion_results.empty:
                    self.logger.info(
                        f"Found conversion results with {len(conversion_results)} rows"
                    )
                    self.logger.info("Generating conversion reports...")
                    report_paths = self.generate_conversion_reports(conversion_results)
                    self.logger.info(f"Generated {len(report_paths)} report files")
                else:
                    self.logger.warning("No conversion results found")
            else:
                self.logger.warning("No response factors found")
        else:
            self.logger.info("Conversion analysis is disabled in config")

        # 3. Generate visualizations
        self.generate_visualizations()

        # 4. Extract sample data
        samples = self.extract_sample_data()

        # 5. Generate compound reports
        report_paths = self.generate_compound_reports(samples)

        self.logger.info("\nAnalysis and visualization complete!")
        self.logger.info(f"Reports saved to: {self.report_dir}")

        return report_paths

    def get_compound_count(self, row, compound_type):
        """
        Get the number of compounds of a specific type by scanning column headers

        Args:
            row: DataFrame row containing sample data
            compound_type: Type of compound (e.g., "reactant", "product")

        Returns:
            Count of compounds detected
        """
        prefix = f"{compound_type}-"
        max_index = 0

        for col in row.index:
            if col.startswith(prefix):
                # Extract index after the prefix (e.g., "reactant-3" -> 3)
                try:
                    parts = col[len(prefix) :].split("-")
                    if parts[0].isdigit():
                        index = int(parts[0])
                        max_index = max(max_index, index)
                except (ValueError, IndexError):
                    continue

        return max_index

    def desalt_smiles(self, smiles: str, return_details: bool = False):
        """Strips salts from a SMILES string by returning the largest molecular fragment."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        try:
            # Skip empty or non-string inputs
            if not isinstance(smiles, str) or not smiles.strip():
                return (smiles, False, []) if return_details else smiles

            # Quick check - if no periods, not a salt
            if "." not in smiles:
                return (smiles, False, []) if return_details else smiles

            # Convert to RDKit molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.logger.warning(f"Could not parse SMILES: {smiles}")
                return (smiles, False, []) if return_details else smiles

            # Calculate molecular weight of original compound
            original_mw = Descriptors.MolWt(mol)
            self.logger.info(
                f"Original compound MW: {original_mw:.2f} Da, SMILES: {smiles}"
            )

            # Get fragments
            fragments = Chem.GetMolFrags(mol, asMols=True)
            if len(fragments) <= 1:
                # No salts present
                return (smiles, False, []) if return_details else smiles

            # Find the largest fragment by molecular weight
            fragment_weights = [Descriptors.MolWt(frag) for frag in fragments]

            # Log all fragments
            for i, (frag, weight) in enumerate(zip(fragments, fragment_weights)):
                frag_smiles = Chem.MolToSmiles(frag)
                self.logger.info(
                    f"  Fragment {i+1}: MW {weight:.2f} Da, SMILES: {frag_smiles}"
                )

            largest_idx = fragment_weights.index(max(fragment_weights))
            main_fragment = fragments[largest_idx]

            # Convert to canonical SMILES
            desalted_smiles = Chem.MolToSmiles(main_fragment)

            self.logger.info(
                f"Selected largest fragment: MW {fragment_weights[largest_idx]:.2f} Da"
            )
            self.logger.info(f"Removed salts from {smiles} -> {desalted_smiles}")

            # Get salt fragments
            salt_fragments = []
            for i, frag in enumerate(fragments):
                if i != largest_idx:
                    salt_smiles = Chem.MolToSmiles(frag)
                    salt_fragments.append(salt_smiles)

            if salt_fragments:
                self.logger.info(f"Salt fragments: {salt_fragments}")

            if return_details:
                return desalted_smiles, True, salt_fragments
            else:
                return desalted_smiles

        except Exception as e:
            self.logger.error(f"Error stripping salts from {smiles}: {str(e)}")
            return (smiles, False, []) if return_details else smiles

    def calculate_response_factors(self):
        """
        Calculate response factors between reference starting materials and internal standards

        The response factor is defined as:
        RF = (reactant_signal/reactant_conc) / (is_signal/is_conc)

        Returns:
            Dictionary mapping compound pairs to their response factors
        """
        self.logger.info("Calculating response factors from reference compounds...")

        # Check if conversion analysis is enabled
        if "conversion" not in self.config or not self.config["conversion"].get(
            "enabled", False
        ):
            self.logger.info("Conversion analysis not enabled in config")
            return {}

        # Get reference configuration
        ref_config = self.config["conversion"]["reference"]
        ref_csv_path = ref_config["csv_file"]
        ref_data_dir = ref_config["data_dir"]
        ref_default_conc = ref_config.get("default_concentration_uM", 100.0)

        # Get internal standards configuration
        is_config = self.config["conversion"]["internal_standards"]
        is_default_conc = is_config.get("default_concentration_uM", 50.0)
        is_standards = is_config.get("standards", [])

        # Load reference data
        try:
            ref_data = pd.read_csv(ref_csv_path)
            self.logger.info(
                f"Loaded reference data from {ref_csv_path}: {len(ref_data)} samples"
            )
        except Exception as e:
            self.logger.error(f"Error loading reference data: {str(e)}")
            return {}

        # Storage for response factors
        response_factors = {}

        # Process each reference sample
        for idx, row in ref_data.iterrows():
            sample_id = row.get("sample-ID", f"reference_{idx}")
            self.logger.info(f"Processing reference sample {sample_id}")

            # Find mzML file
            mzML_filename = row.get("mzML-filename")
            mzML_filepath = self.find_mzml_file(mzML_filename, ref_data_dir)
            if not mzML_filepath:
                self.logger.warning(
                    f"Could not find mzML file for reference sample {sample_id} and filename {mzML_filename}"
                )
                continue

            # Process reactants in this reference
            reactant_count = self.get_compound_count(row, "reactant")
            for i in range(1, reactant_count + 1):
                reactant_col = f"reactant-{i}"
                reactant_name = row.get(reactant_col, f"Reactant-{i}")
                reactant_smiles = row.get(reactant_col)

                # Skip if no SMILES
                if not isinstance(reactant_smiles, str) or not reactant_smiles.strip():
                    continue

                # Desalt SMILES
                reactant_smiles = self.desalt_smiles(reactant_smiles)

                # Get concentration (from CSV or default)
                conc_col = f"reactant-{i}-concentration-uM"
                if conc_col in row.index and not pd.isna(row[conc_col]):
                    reactant_conc = float(row[conc_col])
                else:
                    reactant_conc = ref_default_conc

                self.logger.info(
                    f"Reference reactant: {reactant_name}, Conc: {reactant_conc} μM"
                )

                # Process each internal standard
                for is_info in is_standards:
                    is_name = is_info["name"]
                    is_smiles = is_info["smiles"]
                    is_modes = is_info.get("modes", ["Positive", "Negative"])
                    is_conc = is_info.get("concentration_uM", is_default_conc)
                    is_custom_mw = is_info.get("molecular_weight", None)

                    self.logger.info(
                        f"Using internal standard from config: {is_name}, Conc: {is_conc} μM"
                    )
                    if is_custom_mw is not None:
                        self.logger.info(
                            f"Using custom molecular weight: {is_custom_mw} Da"
                        )

                    # Just check if this internal standard is actually in the reference sample
                    # (we still need to know which standards are present in each sample)
                    is_present = False
                    for j in range(1, 10):  # Check up to 10 internal standards
                        is_col = f"internal-std-{j}"
                        if is_col in row.index and row.get(is_col) == is_name:
                            is_present = True
                            break

                    if not is_present:
                        self.logger.warning(
                            f"Internal standard {is_name} not found in reference sample {sample_id}"
                        )
                        continue

                    # Calculate response factors for each mode
                    for mode in is_modes:
                        if mode not in self.config["parameters"]["modes"]:
                            continue

                        # Extract signals for reactant and IS
                        # This is not getting the correct EIC signal, need to fix this!
                        reactant_signal = self._extract_compound_signal(
                            mzML_filepath, reactant_smiles, mode
                        )
                        is_signal = self._extract_compound_signal(
                            mzML_filepath, is_smiles, mode, custom_mw=is_custom_mw
                        )

                        # Calculate response factor if we have valid signals
                        if reactant_signal > 0 and is_signal > 0:
                            # RF = (reactant_signal/reactant_conc) / (is_signal/is_conc)
                            # This is the relative response factor between the two compounds
                            rf = (reactant_signal / reactant_conc) / (
                                is_signal / is_conc
                            )

                            # Generate a unique key for this reactant/IS/mode combination
                            rf_key = f"{reactant_name}_{is_name}_{mode}"

                            # Store the response factor with metadata
                            response_factors[rf_key] = {
                                "reactant_name": reactant_name,
                                "reactant_smiles": reactant_smiles,
                                "is_name": is_name,
                                "is_smiles": is_smiles,
                                "mode": mode,
                                "response_factor": rf,
                                "reactant_signal": reactant_signal,
                                "is_signal": is_signal,
                                "reactant_conc": reactant_conc,
                                "is_conc": is_conc,
                                "sample_id": sample_id,
                            }

                            self.logger.info(
                                f"Calculated response factor for {reactant_name}/{is_name} ({mode}): {rf:.4f}"
                            )

        self.logger.info(f"Calculated {len(response_factors)} response factors")
        if not response_factors:
            self.logger.warning(
                "No response factors were calculated - check reference data CSV and mzML files"
            )
        return response_factors

    def _extract_compound_signal(self, mzML_filepath, smiles, mode, custom_mw=None):
        """
        Extract EIC area for a compound from an MS file

        Args:
            mzML_filepath: Path to mzML file
            smiles: SMILES string of compound to analyze
            mode: Ionization mode (Positive or Negative)
            custom_mw: Optional custom molecular weight to use instead of calculated value

        Returns:
            EIC area or 0 if not found
        """
        try:
            # Log compound details
            self.logger.info(
                f"Extracting signal for SMILES: {smiles} in {mode} mode with custom MW: {custom_mw}"
            )

            # Calculate or use provided molecular weight
            if custom_mw is not None:
                mw = custom_mw
                self.logger.info(f"Using custom molecular weight: {mw:.2f} Da")
            else:
                # Calculate from SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    self.logger.info(f"Calculated molecular weight: {mw:.2f} Da")
                else:
                    self.logger.warning(f"Could not calculate MW from SMILES: {smiles}")
                    return 0

            # Calculate m/z values to search for
            if mode == "Positive":
                self.logger.info(f"Expected m/z (M+H)+: {mw + 1.007825:.4f}")
            else:
                self.logger.info(f"Expected m/z (M-H)-: {mw - 1.007825:.4f}")

            # Create analyzer for this mode
            analyzer = AnalyseSpectrum(mzMLfilepath=mzML_filepath, mode=mode)

            # Default ions based on mode - need to fix this! Must come from config or csv
            ions_to_add = ["[H]"]

            # Run analysis with custom MW
            analyzer.analyse(
                compoundsmiles=smiles,
                ionstoadd=ions_to_add,
                tolerance=1,
                custom_mw=custom_mw,  # Pass the custom MW
            )

            # Calculate the EIC area - this properly handles the signal extraction
            eic_area = analyzer.calculate_eic_area()
            self.logger.info(f"Calculated EIC area: {eic_area:.2f}")

            return eic_area

        except Exception as e:
            self.logger.error(
                f"Error extracting signal for {smiles} in {mode} mode: {str(e)}"
            )
            return 0

    def calculate_conversions(self, response_factors):
        """
        Calculate conversions for each sample based on response factors

        Args:
            response_factors: Dictionary of response factors between reactants and internal standards

        Returns:
            DataFrame with conversion results
        """
        self.logger.info("Calculating conversions for samples...")

        conversion_results = []

        # Process each sample first
        for idx, sample_data in self.batch_data.iterrows():
            sample_id = sample_data.get("sample-ID", f"sample_{idx}")
            self.logger.info(f"Processing conversion for sample: {sample_id}")

            # Get mzML file for this sample
            mzML_filename = sample_data.get("mzML-filename")
            mzML_filepath = self.find_mzml_file(mzML_filename)
            if not mzML_filepath:
                self.logger.warning(f"Could not find mzML file for sample {sample_id}")
                continue

            # Find which internal standards are present in this sample
            present_standards = {}
            for j in range(1, 10):  # Check up to 10 internal standards
                is_col = f"internal-std-{j}"
                if is_col in sample_data.index and not pd.isna(sample_data[is_col]):
                    is_name = sample_data[is_col]
                    present_standards[is_name] = (
                        j  # Store position index for potential concentration lookup
                    )
                    self.logger.info(
                        f"Found internal standard: {is_name} in sample {sample_id}"
                    )

            if not present_standards:
                self.logger.warning(
                    f"No internal standards found in sample {sample_id}"
                )
                continue

            # Find available reactants in this sample
            reactant_count = self.get_compound_count(sample_data, "reactant")

            # Process each reactant/internal standard combination that has a response factor
            for is_name, is_idx in present_standards.items():
                # Get internal standard info from config
                is_config = None
                is_custom_mw = None
                is_smiles = None
                for is_info in self.config["conversion"]["internal_standards"].get(
                    "standards", []
                ):
                    if is_info["name"] == is_name:
                        is_config = is_info
                        is_custom_mw = is_info.get("molecular_weight", None)
                        is_smiles = is_info["smiles"]
                        is_conc = is_info.get(
                            "concentration_uM",
                            self.config["conversion"]["internal_standards"].get(
                                "default_concentration_uM", 100.0
                            ),
                        )
                        break

                if is_config is None:
                    self.logger.warning(
                        f"Internal standard {is_name} not found in configuration"
                    )
                    continue

                # Extract internal standard signal from this sample
                is_signal = self._extract_compound_signal(
                    mzML_filepath, is_smiles, "Positive", custom_mw=is_custom_mw
                )
                if is_signal == 0:
                    self.logger.warning(
                        f"No signal detected for internal standard {is_name} in sample {sample_id}"
                    )
                    continue

                # Process each reactant in this sample
                for i in range(1, reactant_count + 1):
                    reactant_col = f"reactant-{i}"
                    reactant_name = sample_data.get(reactant_col, f"Reactant-{i}")
                    reactant_smiles = sample_data.get(reactant_col)

                    # Skip if no SMILES
                    if (
                        not isinstance(reactant_smiles, str)
                        or not reactant_smiles.strip()
                    ):
                        continue

                    # Desalt SMILES
                    reactant_smiles = self.desalt_smiles(reactant_smiles)

                    # Check if we have a response factor for this pair
                    rf_key_pos = f"{reactant_name}_{is_name}_Positive"
                    rf_key_neg = f"{reactant_name}_{is_name}_Negative"

                    for mode, rf_key in [
                        ("Positive", rf_key_pos),
                        ("Negative", rf_key_neg),
                    ]:
                        if rf_key not in response_factors:
                            continue

                        # Get response factor data
                        rf_data = response_factors[rf_key]
                        response_factor = rf_data["response_factor"]

                        # Get reactant signal (if exists in batch data)
                        signal_column = f"{reactant_col}-EIC-area-{mode}"
                        if signal_column in sample_data.index and not pd.isna(
                            sample_data[signal_column]
                        ):
                            reactant_signal = sample_data[signal_column]

                            # Calculate conversion
                            concentration_reactant = (
                                ((reactant_signal / is_signal) * is_conc)
                                / response_factor
                                if is_signal > 0 and response_factor > 0
                                else 0
                            )

                            # Store result
                            conversion_results.append(
                                {
                                    "sample_id": sample_id,
                                    "reactant_name": reactant_name,
                                    "reactant_idx": i,
                                    "internal_standard": is_name,
                                    "mode": mode,
                                    "reactant_signal": reactant_signal,
                                    "is_signal": is_signal,
                                    "response_factor": response_factor,
                                    "is_concentration_uM": is_conc,
                                    "concentration_reactant_uM": concentration_reactant,
                                }
                            )

                            # Also add to batch data
                            col_name = f"reactant-{i}-concentration-{mode}-{is_name}"
                            self.batch_data.loc[idx, col_name] = concentration_reactant

                            self.logger.info(
                                f"Calculated concentration for {reactant_name} using {is_name}: {concentration_reactant:.2f}%"
                            )

        # Create DataFrame from results
        if conversion_results:
            df = pd.DataFrame(conversion_results)
            self.logger.info(
                f"Created conversion results DataFrame with {len(df)} rows"
            )
            return df
        else:
            self.logger.warning("No conversion results were calculated")
            return pd.DataFrame()

    # Helper function to extract reactant index from column name
    def _get_reactant_idx(self, column_name):
        """Extract reactant index from column name like 'reactant-1-max-EIC-signal-Positive'"""
        match = re.search(r"reactant-(\d+)", column_name)
        if match:
            return match.group(1)
        return None


# Configure logging
logs_dir = "/Users/bvh64415/myrepos/mscheck/logs"

# Create timestamped log filename
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = os.path.join(logs_dir, f"mscheck_{timestamp}.log")

setup_logger(name="MSCheck", level="INFO", log_file=log_file, use_colors=True)
logger = get_logger("MSCheck")  # Use get_logger for consistency with your module
logger.info(f"Logging to file: {log_file}")

# Initialize with config
analyzer = BulkAnalyser(
    "/Users/bvh64415/myrepos/mscheck/tests/testdata/bulk-test/mscheck_config_no_conversion.yaml"
)

# # Run specific steps as needed
# analyzer.load_data()
# analyzer.process_samples(["reactant", "product"], ["Positive"])
# response_factors = analyzer.calculate_response_factors()
# conversion_results = analyzer.calculate_conversions(response_factors)

# Or run the complete workflow
report_paths = analyzer.run_complete_workflow()

logger.info(f"Analysis complete! Generated {len(report_paths)} reports")
