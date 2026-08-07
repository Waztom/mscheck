from datetime import datetime
import os
import gc
import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors
import yaml

# Local imports
from .analyseMS import AnalyseMS
from .analyseUV import AnalyseUV
from .report import MSReport
from .heatmap import MSHeatmapGenerator
from .utils import monitor_memory, desalt_smiles, standardise_compound


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

    
    def analyse_batch(self, batch_size=20, analysis_types=None, modes=None, tolerance=None):
        """
        Process samples in smaller batches to reduce memory usage
        
        Args:
            batch_size: Number of samples to process in each batch
            analysis_types: List of analysis types (default: from config)
            modes: List of ionization modes (default: from config)
            tolerance: Mass tolerance in ppm (default: from config)
            
        Returns:
            DataFrame with analysis results
        """
        # Use config parameters if not specified
        if analysis_types is None:
            analysis_types = self.config["parameters"]["analysis_types"]
        if modes is None:
            modes = self.config["parameters"]["modes"]
        if tolerance is None:
            tolerance = self.config["parameters"]["tolerance"]
            
        self.logger.info(f"Processing samples in batches of {batch_size}...")
        monitor_memory("Before batch processing", self.logger)
        
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
        
        # Calculate number of batches
        total_samples = len(self.batch_data)
        num_batches = (total_samples + batch_size - 1) // batch_size  # Ceiling division
        
        self.logger.info(f"Processing {total_samples} samples in {num_batches} batches")
        
        # Process each batch
        for batch_idx in range(num_batches):
            batch_start = batch_idx * batch_size
            batch_end = min(batch_start + batch_size, total_samples)
            
            self.logger.info(f"Processing batch {batch_idx+1}/{num_batches} (samples {batch_start+1}-{batch_end})")
            
            # Process just this subset of samples
            batch_indices = list(range(batch_start, batch_end))
            
            # Process each row in the batch
            for batch_index in batch_indices:
                row = self.batch_data.iloc[batch_index]
                try:
                    sample_id = str(row.get("sample-ID", f"sample_{batch_index}"))
                    mzML_filename = row["mzML-filename"]
                    self.logger.info(f"Processing sample {sample_id} ({batch_index + 1}/{total_samples})")

                    # Find mzML file
                    mzML_filepath = self.find_mzml_file(mzML_filename)
                    if mzML_filepath is None:
                        continue

                    # Precompute per-analysis-type row info (independent of mode) so we can
                    # reuse a single AnalyseMS instance per (sample, mode) across every
                    # analysis_type / compound instead of re-parsing the mzML each time.
                    type_infos = []
                    for analysis_type in analysis_types:
                        no_compounds = self.get_compound_count(row, analysis_type)
                        if no_compounds <= 0:
                            continue
                        ions_add_column = f"{analysis_type}-ions-to-add"
                        ions_sub_column = f"{analysis_type}-ions-to-sub"
                        tolerance_column = f"{analysis_type}-match-tolerance"
                        type_infos.append({
                            "type": analysis_type,
                            "no_compounds": no_compounds,
                            "ions_to_add": self.parse_ions(row.get(ions_add_column, "")) or ["[H]"],
                            "ions_to_sub": self.parse_ions(row.get(ions_sub_column, "")),
                            "match_tolerance": int(row.get(tolerance_column, tolerance)),
                        })

                    # One AnalyseMS per (sample, mode); reused across every compound.
                    for mode in modes:
                        analysis_obj = None
                        try:
                            analysis_obj = AnalyseMS(mzMLfilepath=mzML_filepath, mode=mode)
                            report_subdir = os.path.join(self.report_dir, mode)

                            for info in type_infos:
                                analysis_type = info["type"]
                                for compound_idx in range(info["no_compounds"]):
                                    compound_column = f"{analysis_type}-{compound_idx + 1}"
                                    if compound_column not in row:
                                        continue

                                    analysis_type_smiles = row[compound_column]

                                    # Special handling for internal standards
                                    if analysis_type == "internal-std":
                                        std_name = analysis_type_smiles
                                        for is_info in self.config.get("conversion", {}).get("internal_standards", {}).get("standards", []):
                                            if is_info.get("name") == std_name:
                                                analysis_type_smiles = is_info.get("smiles")
                                                break

                                    if not isinstance(analysis_type_smiles, str) or not analysis_type_smiles.strip():
                                        continue

                                    analysis_type_smiles = desalt_smiles(analysis_type_smiles)

                                    try:
                                        total_analyses_attempted += 1

                                        analysis_obj.analyse(
                                            compoundsmiles=analysis_type_smiles,
                                            ionstoadd=info["ions_to_add"],
                                            ionstosub=info["ions_to_sub"],
                                            tolerance=info["match_tolerance"],
                                        )

                                        # Calculate EIC area
                                        eic_area = analysis_obj.calculate_eic_area()

                                        # Get max RT for an ion match
                                        if "max_rt" in analysis_obj.analysedata and analysis_obj.analysedata["max_rt"] is not None:
                                            rt_column = f"{analysis_type}-{compound_idx + 1}-RT-{mode}"
                                            rt_value = analysis_obj.analysedata["max_rt"]
                                            self.batch_data.loc[batch_index, rt_column] = rt_value
                                            self.logger.debug(f"Stored max RT {rt_value:.2f} min for {analysis_type}-{compound_idx + 1} in {mode} mode")

                                        # Create report
                                        report_name = f"{sample_id}_{analysis_type}_{compound_idx + 1}"
                                        analysis_obj.create_report(folder=report_subdir, compound_name=report_name)

                                        # Update batch data
                                        signal_column = f"{analysis_type}-{compound_idx + 1}-max-EIC-signal-{mode}"
                                        mz_match_column = f"{analysis_type}-{compound_idx + 1}-max-mz-match-{mode}"
                                        ions_column = f"{analysis_type}-{compound_idx + 1}-ions-matched-{mode}"
                                        eic_area_column = f"{analysis_type}-{compound_idx + 1}-EIC-area-{mode}"

                                        self.batch_data.loc[batch_index, signal_column] = analysis_obj.analysedata["max_EIC_signal"]
                                        self.batch_data.loc[batch_index, mz_match_column] = str(analysis_obj.analysedata["max_mz_match"])
                                        self.batch_data.loc[batch_index, ions_column] = str(analysis_obj.analysedata["ions"])
                                        self.batch_data.loc[batch_index, eic_area_column] = eic_area

                                        successful_analyses += 1

                                    except Exception as e:
                                        self.logger.error(f"Error analyzing {analysis_type}-{compound_idx + 1} in {mode}: {str(e)}")
                        finally:
                            if analysis_obj is not None:
                                del analysis_obj

                    self.processed_samples += 1

                    # Process UV data for this sample if MS analysis was successful
                    try:
                        self.logger.info(f"Processing UV data for sample {sample_id}")
                        
                        # Find mzML file path (we already have it from MS analysis)
                        if mzML_filepath:
                            # Create UV analyzer
                            uv_analyzer = None
                            try:
                                uv_analyzer = AnalyseUV(mzML_filepath)
                                
                                # Skip if no UV data available
                                if not uv_analyzer.has_uv:
                                    self.logger.warning(f"No UV data available for sample {sample_id}")
                                else:
                                    self.logger.info(f"Found UV data with {len(uv_analyzer.chromatograms)} chromatograms")
                                    
                                    # Process UV data for each compound that had successful MS analysis
                                    for analysis_type in analysis_types:
                                        # Get number of compounds of this type
                                        no_compounds = self.get_compound_count(row, analysis_type)
                                        if no_compounds <= 0:
                                            continue
                                        
                                        for compound_idx in range(no_compounds):
                                            compound_col = f"{analysis_type}-{compound_idx + 1}"
                                            compound_name = f"{analysis_type}-{compound_idx + 1}"
                                            
                                            # Skip if compound not in data
                                            if compound_col not in row:
                                                continue
                                            
                                            # Get retention time from MS analysis
                                            rt_from_ms = None
                                            for mode in modes:
                                                rt_col = f"{analysis_type}-{compound_idx + 1}-RT-{mode}"
                                                if rt_col in self.batch_data.columns and not pd.isna(self.batch_data.loc[batch_index, rt_col]):
                                                    rt_from_ms = float(self.batch_data.loc[batch_index, rt_col])
                                                    self.logger.info(f"Using RT {rt_from_ms:.2f} min from {mode} mode for UV analysis")
                                                    break
                                                    
                                            # If RT was found, perform UV analysis at that retention time
                                            if rt_from_ms is not None:
                                                self.logger.info(f"Analyzing UV data for {compound_name} at RT {rt_from_ms:.2f} min")
                                                
                                                # Set tolerance based on peak width
                                                tolerance = 0.2  # Default 0.2 min window
                                                
                                                # Analyze UV data at the retention time
                                                results = uv_analyzer.analyse(retention_time=rt_from_ms, tolerance=tolerance)
                                            
                                                # Store results if peak found
                                                if results['found']:
                                                    # Update batch data
                                                    max_area_col = f"{analysis_type}-{compound_idx + 1}-UV-max-area"
                                                    optimal_wl_col = f"{analysis_type}-{compound_idx + 1}-UV-optimal-wavelength"
                                                    max_intensity_col = f"{analysis_type}-{compound_idx + 1}-UV-max-intensity"
                                                    
                                                    self.batch_data.loc[batch_index, max_area_col] = results['max_area']
                                                    self.batch_data.loc[batch_index, optimal_wl_col] = results['optimal_wavelength']
                                                    self.batch_data.loc[batch_index, max_intensity_col] = results['max_intensity']
                                                    
                                                    self.logger.info(f"UV analysis for {compound_name}: wavelength={results['optimal_wavelength']} nm, Area={results['max_area']:.2f}")
                                                else:
                                                    self.logger.warning(f"No UV peak found for {compound_name} at RT {rt_from_ms:.2f} min")
                                            else:
                                                self.logger.warning(f"No retention time found for {compound_name}")
                            
                            finally:
                                # Clean up analyzer to free memory
                                if uv_analyzer is not None:
                                    del uv_analyzer
                                    gc.collect()
                                    
                    except Exception as e:
                        self.logger.error(f"Error analyzing UV data for sample {sample_id}: {str(e)}")

                except Exception as e:
                    self.logger.error(f"Error processing sample at index {batch_index}: {str(e)}")

                # Reclaim memory once per sample, after MS + UV analyzers are gone.
                gc.collect()

            # Save intermediate results after each batch
            intermediate_path = os.path.join(self.report_dir, f"intermediate_batch_{batch_idx+1}.csv")
            self.batch_data.to_csv(intermediate_path, index=False)
            self.logger.info(f"Saved intermediate results for batch {batch_idx+1} to {intermediate_path}")

            gc.collect()
            monitor_memory(f"After processing batch {batch_idx+1}/{num_batches}", self.logger)
        
        self.logger.info(f"Batch processing complete: {successful_analyses}/{total_analyses_attempted} analyses successful")
        self.logger.info(f"Processed {self.processed_samples}/{total_samples} samples")
        
        # One final memory cleanup
        gc.collect()
        monitor_memory("After batch processing complete", self.logger)
        
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

        return [s.strip() for s in ion_str.split(",") if s.strip()]

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
        analysis_types=None,
        modes=None,
        tolerance=None,
        batch_size=20
    ):
        """
        Process all samples and generate individual reports using batch processing
        
        Args:
            analysis_types: List of compound types to analyze (e.g., product, reactant)
            modes: List of ionization modes (e.g., Positive, Negative)
            tolerance: Default mass tolerance in ppm
            batch_size: Number of samples to process in each batch
            
        Returns:
            DataFrame with analysis results
        """
        if self.batch_data is None:
            self.load_data()

        # Process the samples in batches with the specified parameters
        self.analyse_batch(
            batch_size=batch_size,
            analysis_types=analysis_types,
            modes=modes,
            tolerance=tolerance
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

        temp_analyzer = None
        try:
            # Find mzML file
            mzML_filename = row["mzML-filename"]
            mzML_filepath = self.find_mzml_file(mzML_filename)
            if mzML_filepath is not None:
                # Create analyzer to get RT/TIC
                mode = "Positive"
                temp_analyzer = AnalyseMS(mzMLfilepath=mzML_filepath, mode=mode)
                
                # Create copies instead of references to ensure data persists after analyzer is deleted
                if "RT" in temp_analyzer.MSdata:
                    samples[sample_id]["RT_values"] = temp_analyzer.MSdata["RT"].copy() 
                if "TIC" in temp_analyzer.MSdata:
                    samples[sample_id]["TIC_values"] = temp_analyzer.MSdata["TIC"].copy()

                # Store MZ data if available - make a copy to avoid reference issues
                if "mz_data" in temp_analyzer.MSdata:
                    samples[sample_id]["all_mz_data"] = temp_analyzer.MSdata["mz_data"].copy()
                    self.logger.info(f"Added complete MZ data to sample {sample_id}")
        except Exception as e:
            self.logger.error(f"Error loading RT/TIC for sample {sample_id}: {str(e)}")
        finally:
            # Clean up analyzer object
            if temp_analyzer is not None:
                del temp_analyzer
                import gc
                gc.collect()

    def _extract_compound_data(self, samples, sample_id, row, compound_type):
        """Extract compound data for a specific compound type"""
        # Get number of compounds by scanning column headers
        no_compounds = self.get_compound_count(row, compound_type)

        # Process each compound
        for compound_idx in range(no_compounds):
            compound_results = self._extract_single_compound(
                sample_id, row, compound_type, compound_idx
            )

            # Handle the case where we get a list of results (one per mode)
            if isinstance(compound_results, list):
                for result in compound_results:
                    if result:  # Skip None results
                        samples[sample_id]["compounds"].append(result)
            elif compound_results:  # For backward compatibility
                samples[sample_id]["compounds"].append(compound_results)

    def _extract_single_compound(self, sample_id, row, compound_type, compound_idx):
        """Extract data for a single compound in all configured modes"""
        # Get compound information
        compound_column = f"{compound_type}-{compound_idx + 1}"
        smiles_value = row.get(compound_column)

        # Skip if no SMILES or name
        if not isinstance(smiles_value, str) or not smiles_value.strip():
            return None
            
        # Special handling for internal standards - convert name to SMILES
        if compound_type == "internal-std":
            std_name = smiles_value  # This is actually the name, not SMILES
            smiles_found = False
            
            # Look up SMILES from config
            for is_info in (self.config.get("conversion", {})
                            .get("internal_standards", {})
                            .get("standards", [])):
                if is_info.get("name") == std_name:
                    smiles_value = is_info.get("smiles")
                    smiles_found = True
                    self.logger.info(f"Found SMILES for internal standard {std_name}: {smiles_value}")
                    break
                    
            if not smiles_found:
                self.logger.warning(f"No SMILES found in config for internal standard: {std_name}")
                return None

        # Desalt SMILES
        smiles_value = desalt_smiles(smiles_value)
        
        # Create RDKit molecule
        try:
            mol = Chem.MolFromSmiles(smiles_value)
            if mol is None:
                self.logger.warning(f"Could not create molecule from SMILES: {smiles_value}")
                return None
        except Exception as e:
            self.logger.error(f"Error creating molecule: {str(e)}")
            return None

        # Get available modes from config
        modes = self.config.get("parameters", {}).get("modes", ["Positive"])
        
        # Process each mode and collect all results
        compound_results = []
        
        for mode in modes:
            try:
                # Construct column names based on mode
                signal_column = f"{compound_type}-{compound_idx + 1}-max-EIC-signal-{mode}"
                eic_area_column = f"{compound_type}-{compound_idx + 1}-EIC-area-{mode}" 
                mz_match_column = f"{compound_type}-{compound_idx + 1}-max-mz-match-{mode}"
                
                # Get MS data for this mode
                rt_max, intensity_max, eic_data, mz_strongest, all_mz_data = (
                    self._get_compound_ms_data(
                        sample_id, row, compound_type, compound_idx, smiles_value, mode
                    )
                )
                
                # Create result for this mode
                mode_result = {
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
                
                # Add this mode's result to the collection
                compound_results.append(mode_result)
                
            except Exception as e:
                self.logger.error(f"Error processing {compound_type} {compound_idx+1} in {mode} mode: {str(e)}")
    
        # Also need to update _extract_compound_data to handle lists:
        return compound_results

    def _get_compound_ms_data(
        self, sample_id, row, compound_type, compound_idx, smiles_value, mode
    ):
        """Get mass spectrum data for a compound with proper analyzer cleanup"""
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

        temp_analyzer = None
        try:
            # Get mzML file
            mzML_filename = row["mzML-filename"]
            mzML_filepath = self.find_mzml_file(mzML_filename)

            # Create analyzer and run analysis
            temp_analyzer = AnalyseMS(mzMLfilepath=mzML_filepath, mode=mode)

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

            # Extract RT max and intensity max - make copies of arrays
            if ("RT" in temp_analyzer.analysedata 
                    and len(temp_analyzer.analysedata["RT"]) > 0):
                rt_values = temp_analyzer.analysedata["RT"][0]
                tic_values = temp_analyzer.analysedata["TIC"][0]
                
                # Make copies if they're numpy arrays
                if isinstance(rt_values, np.ndarray):
                    rt_values = rt_values.copy()
                if isinstance(tic_values, np.ndarray):
                    tic_values = tic_values.copy()
                    
                if len(rt_values) > 0 and len(tic_values) > 0:
                    max_idx = np.argmax(tic_values)
                    rt_max = rt_values[max_idx]
                    intensity_max = tic_values[max_idx]

            # Extract EIC data - make a copy
            if ("EIC_data" in temp_analyzer.analysedata 
                    and len(temp_analyzer.analysedata["EIC_data"]) > 0):
                eic_data = temp_analyzer.analysedata["EIC_data"][0]
                if isinstance(eic_data, np.ndarray):
                    eic_data = eic_data.copy()

            # Extract mass spectrum data - make copies
            if ("mz_strongest" in temp_analyzer.analysedata 
                    and len(temp_analyzer.analysedata["mz_strongest"]) > 0):
                mz_masses, mz_intensities, _ = temp_analyzer.analysedata["mz_strongest"][0]
                
                # Make copies if they're numpy arrays
                if isinstance(mz_masses, np.ndarray):
                    mz_masses = mz_masses.copy()
                if isinstance(mz_intensities, np.ndarray):
                    mz_intensities = mz_intensities.copy()
                    
                mz_strongest = (mz_masses, mz_intensities)

            # Get full MZ data - handle different types
            if hasattr(temp_analyzer, "MSdata") and "mz_data" in temp_analyzer.MSdata:
                mz_data = temp_analyzer.MSdata["mz_data"]
                
                # Handle different types of mz_data
                if isinstance(mz_data, dict):
                    # If it's a dictionary, copy each item
                    all_mz_data = {}
                    for key, value in mz_data.items():
                        if isinstance(value, np.ndarray):
                            all_mz_data[key] = value.copy()
                        else:
                            all_mz_data[key] = value
                elif isinstance(mz_data, list):
                    # If it's a list, make a copy of the list
                    all_mz_data = []
                    for item in mz_data:
                        if isinstance(item, np.ndarray):
                            all_mz_data.append(item.copy())
                        else:
                            all_mz_data.append(item)
                else:
                    # For other types, just store directly
                    all_mz_data = mz_data

        except Exception as e:
            self.logger.error(f"Error extracting MS data: {str(e)}")
        
        finally:
            # Always clean up the analyzer object, even if an error occurs
            if temp_analyzer is not None:
                del temp_analyzer
                import gc
                gc.collect()

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
        Generate reports and visualizations for conversion results with improved statistics
        
        Args:
            conversion_df: DataFrame with conversion results
            
        Returns:
            Dictionary of report paths
        """
        self.logger.info("Generating conversion reports with improved statistics...")
        
        if conversion_df is None or conversion_df.empty:
            self.logger.warning("No conversion data available for reports")
            return {}
            
        # Create report directory
        report_dir = os.path.join(self.report_dir, "conversion_reports")
        os.makedirs(report_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_paths = {}
        
        # Make a copy to avoid modifying the original
        working_df = conversion_df.copy()
        
        # Add UV area data to conversion summary
        self.logger.info("Adding UV area data to conversion summary...")

        # Extract UV area data from batch_data
        uv_area_columns = [col for col in self.batch_data.columns if '-UV-max-area' in col]
        if uv_area_columns:
            self.logger.info(f"Found {len(uv_area_columns)} UV area columns to include in reports")
            
            # Create mapping between sample IDs in working_df and indices in batch_data
            sample_id_to_index = {}
            for idx, row in self.batch_data.iterrows():
                sample_id = row.get("sample-ID", f"sample_{idx}")
                sample_id_to_index[sample_id] = idx
            
            # For each sample in conversion results, add UV area data
            for idx, row in working_df.iterrows():
                sample_id = row["sample_id"]
                if sample_id in sample_id_to_index:
                    batch_idx = sample_id_to_index[sample_id]
                    
                    # Add each UV area column with a cleaner name
                    for col in uv_area_columns:
                        # Extract compound type and number from column name
                        # Example: "reactant-1-UV-max-area" -> "reactant-1-UVArea"
                        parts = col.split('-UV-')
                        compound_id = parts[0]  # e.g., "reactant-1"
                        
                        # Create new column name
                        new_col = f"{compound_id}-UVArea"
                        
                        # Get the UV area value if it exists
                        if not pd.isna(self.batch_data.loc[batch_idx, col]):
                            working_df.loc[idx, new_col] = self.batch_data.loc[batch_idx, col]
                            
            self.logger.info("UV area data added to conversion summary")
        else:
            self.logger.info("No UV area data found in batch results")

        # Find input concentrations if available by looking at reactant configuration
        input_concentrations = {}
        if "conversion" in self.config and "reactants" in self.config["conversion"]:
            for reactant in self.config["conversion"]["reactants"]:
                if "name" in reactant and "concentration_uM" in reactant:
                    input_concentrations[reactant["name"]] = reactant["concentration_uM"]
        
        # Filter out unrealistic values (concentration ≥ input concentration)
        filtered_records = 0
        rows_to_drop = []

        for idx, row in working_df.iterrows():  # Iterate over working_df directly
            reactant = row["reactant_name"]
            
            # Skip rows with NaN concentration
            if pd.isna(row["concentration_reactant_uM"]):
                continue
                
            concentration = float(row["concentration_reactant_uM"])
            
            # If we know the input concentration, filter out if concentration is higher
            if reactant in input_concentrations:
                input_conc = input_concentrations[reactant]
                if concentration >= input_conc:
                    rows_to_drop.append(idx)
                    filtered_records += 1
                    
        # Drop all identified rows at once
        if rows_to_drop:
            working_df = working_df.drop(rows_to_drop)
            
        if filtered_records > 0:
            self.logger.info(f"Filtered out {filtered_records} records with unrealistic concentrations")
        
        # ===== 1. Generate detailed CSV summary =====
        summary_path = os.path.join(report_dir, f"conversion_summary_{timestamp}.csv")
        working_df.to_csv(summary_path, index=False)
        report_paths["detailed_summary"] = summary_path
        self.logger.info(f"Saved detailed conversion summary to: {summary_path}")
        
        # Define custom aggregation function to handle single-value case
        def calc_concentration(x):
            if len(x) == 1:  # Single value
                return x.iloc[0]
            else:  # Multiple values
                return x.mean()
        
        # ===== 2. Generate summary by sample =====
        try:
            # Group by sample and calculate statistics with improved method
            sample_summary = working_df.groupby("sample_id").agg({
                "concentration_reactant_uM": [
                    calc_concentration,  # Custom function 
                    "std", 
                    "min", 
                    "max",
                    "count"
                ]
            }).reset_index()
            
            # Flatten column names
            sample_summary.columns = [
                "sample_id", 
                "calculated_concentration",  # New name for our custom function
                "std_concentration", 
                "min_concentration", 
                "max_concentration",
                "measurement_count"
            ]

            # Add UV area columns if available
            if uv_area_columns:
                # For each sample in sample_summary, add average UV areas
                for idx, row in sample_summary.iterrows():
                    sample_id = row["sample_id"]
                    if sample_id in sample_id_to_index:
                        batch_idx = sample_id_to_index[sample_id]
                        
                        # Add each UV area column with a cleaner name
                        for col in uv_area_columns:
                            parts = col.split('-UV-')
                            compound_id = parts[0]  # e.g., "reactant-1"
                            new_col = f"{compound_id}-UVArea"
                            
                            if not pd.isna(self.batch_data.loc[batch_idx, col]):
                                sample_summary.loc[idx, new_col] = self.batch_data.loc[batch_idx, col]
                
                self.logger.info("UV area data added to sample summary")  

            sample_path = os.path.join(report_dir, f"sample_summary_{timestamp}.csv")
            sample_summary.to_csv(sample_path, index=False)
            report_paths["sample_summary"] = sample_path
            self.logger.info(f"Saved sample conversion summary to: {sample_path}")
        except Exception as e:
            self.logger.error(f"Error generating sample summary: {str(e)}")
        
        # ===== 3. Generate summary by reactant =====
        try:
            # Group by reactant with improved method
            reactant_summary = working_df.groupby("reactant_name").agg({
                "concentration_reactant_uM": [
                    calc_concentration,  # Custom function
                    "std", 
                    "min", 
                    "max", 
                    "count"
                ]
            }).reset_index()
            
            # Flatten column names
            reactant_summary.columns = [
                "reactant_name",
                "calculated_concentration",  # New name for our custom function
                "std_concentration", 
                "min_concentration", 
                "max_concentration",
                "measurement_count"
            ]
            
            reactant_path = os.path.join(report_dir, f"reactant_summary_{timestamp}.csv")
            reactant_summary.to_csv(reactant_path, index=False)
            report_paths["reactant_summary"] = reactant_path
            self.logger.info(f"Saved reactant conversion summary to: {reactant_path}")
        except Exception as e:
            self.logger.error(f"Error generating reactant summary: {str(e)}")
        
        # ===== 4. Generate summary by internal standard and reactant =====
        try:
            # Group by reactant and internal standard with improved method
            is_summary = working_df.groupby(["reactant_name", "internal_standard", "mode"]).agg({
                "concentration_reactant_uM": [
                    calc_concentration,  # Custom function
                    "std", 
                    "min", 
                    "max", 
                    "count"
                ],
                "response_factor": ["mean"]
            }).reset_index()
            
            # Flatten column names
            is_summary.columns = [
                "reactant_name", 
                "internal_standard", 
                "mode",
                "calculated_concentration",  # New name for our custom function
                "std_concentration", 
                "min_concentration", 
                "max_concentration",
                "measurement_count",
                "mean_response_factor"
            ]
            
            is_path = os.path.join(report_dir, f"is_comparison_{timestamp}.csv")
            is_summary.to_csv(is_path, index=False)
            report_paths["is_summary"] = is_path
            self.logger.info(f"Saved internal standard comparison to: {is_path}")
        except Exception as e:
            self.logger.error(f"Error generating internal standard comparison: {str(e)}")
    
        return report_paths

    def run_complete_workflow(self, batch_size=20):
        """
        Run the complete analysis workflow using parameters from config
        
        Args:
            batch_size: Number of samples to process in each batch
            
        Returns:
            Dictionary mapping sample IDs to their report paths
        """
        # Log initial memory usage
        monitor_memory("Workflow start", self.logger)
        
        # Get parameters from config
        analysis_types = self.config["parameters"]["analysis_types"]
        modes = self.config["parameters"]["modes"]
        tolerance = self.config["parameters"]["tolerance"]

        self.logger.info("Starting complete analysis workflow with batch processing...")
        self.logger.info(f"Analysis types: {', '.join(analysis_types)}")
        self.logger.info(f"Modes: {', '.join(modes)}")
        self.logger.info(f"Default tolerance: {tolerance}")
        self.logger.info(f"Batch size: {batch_size} samples")

        # 1. Process all samples with the specified parameters using batches
        self.process_samples(
            analysis_types=analysis_types, 
            modes=modes, 
            tolerance=tolerance,
            batch_size=batch_size
        )
        monitor_memory("After sample processing", self.logger)

        # 2. Calculate response factors and conversions if enabled
        if "conversion" in self.config and self.config["conversion"].get("enabled", True):
            self.logger.info("Conversion analysis is enabled in config")
            self.logger.info("Calculating response factors for conversion analysis...")

            response_factors = self.calculate_response_factors()
            self.logger.info(f"Found {len(response_factors)} response factors")
            monitor_memory("After response factor calculation", self.logger)

            if response_factors:
                self.logger.info("Calculating conversions using response factors...")
                conversion_results = self.calculate_conversions(response_factors)
                monitor_memory("After conversion calculation", self.logger)

                if not conversion_results.empty:
                    self.logger.info(f"Found conversion results with {len(conversion_results)} rows")
                    self.logger.info("Generating conversion reports...")
                    report_paths = self.generate_conversion_reports(conversion_results)
                    self.logger.info(f"Generated {len(report_paths)} report files")
                    
                    # Clean up large dataframe
                    del conversion_results
                    gc.collect()
                    monitor_memory("After conversion reports", self.logger)
                else:
                    self.logger.warning("No conversion results found")
            else:
                self.logger.warning("No response factors found")
        else:
            self.logger.info("Conversion analysis is disabled in config")

        # 3. Generate visualizations
        self.generate_visualizations()
        monitor_memory("After visualizations", self.logger)

        # 4. Extract sample data
        samples = self.extract_sample_data()
        monitor_memory("After sample data extraction", self.logger)

        # 5. Generate compound reports
        report_paths = self.generate_compound_reports(samples)
        monitor_memory("After compound reports", self.logger)

        self.logger.info("\nAnalysis and visualization complete!")
        self.logger.info(f"Reports saved to: {self.report_dir}")
        
        # Final memory usage
        monitor_memory("Workflow complete", self.logger)

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

                # Standardize the reactant SMILES
                original_smiles = reactant_smiles
                reactant_smiles = desalt_smiles(reactant_smiles)
                reactant_inchi = standardise_compound(reactant_smiles, logger=self.logger)
                self.logger.info(f"Standardized reference reactant: {original_smiles} → {reactant_inchi or 'FAILED'}")

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

                            # When creating the response factor key, use InChI instead of the reactant name:
                           
                            rf_key = f"{reactant_inchi}_{is_name}_{mode}"

                            # Store the response factor with metadata
                            response_factors[rf_key] = {
                                "reactant_name": reactant_name,
                                "reactant_smiles": reactant_smiles,
                                "reactant_inchi": reactant_inchi,
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
        analyzer = None
        try:
            # Log compound details (keep existing logging)
            self.logger.info(f"Extracting signal for SMILES: {smiles} in {mode} mode with custom MW: {custom_mw}")

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
                    # Clean up the molecule object
                    del mol
                else:
                    self.logger.warning(f"Could not calculate MW from SMILES: {smiles}")
                    return 0

            # Keep existing m/z logging
            if mode == "Positive":
                self.logger.info(f"Expected m/z (M+H)+: {mw + 1.007825:.4f}")
            else:
                self.logger.info(f"Expected m/z (M-H)-: {mw - 1.007825:.4f}")

            # Create analyzer for this mode
            analyzer = AnalyseMS(mzMLfilepath=mzML_filepath, mode=mode)

            # Default ions based on mode - keep existing code
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
            self.logger.error(f"Error extracting signal for {smiles} in {mode} mode: {str(e)}")
            return 0
            
        finally:
            # Explicit cleanup to free memory
            if analyzer is not None:
                del analyzer
                import gc
                gc.collect()

    def calculate_conversions(self, response_factors):
        """
        Calculate conversions for each sample with comprehensive reporting of all samples
        
        Args:
            response_factors: Dictionary of response factors from calculate_response_factors
            
        Returns:
            DataFrame with conversion results for all samples, including status information
        """
        self.logger.info("Calculating conversions for samples with comprehensive reporting...")
        self.logger.info(f"Working with {len(response_factors)} response factors")
        
        # Prepare results list that will include ALL samples
        conversion_results = []
        
        # Process each sample in the batch data
        for idx, sample_data in self.batch_data.iterrows():
            sample_id = sample_data.get("sample-ID", f"sample_{idx}")
            self.logger.info(f"Processing sample: {sample_id}")
            
            # Get mzML file path - first potential failure point
            mzML_filename = sample_data.get("mzML-filename")
            mzML_filepath = self.find_mzml_file(mzML_filename)
            
            if not mzML_filepath:
                # Add entry for missing mzML file
                conversion_results.append({
                    "sample_id": sample_id,
                    "reactant_name": "UNKNOWN",
                    "reactant_idx": 0,
                    "internal_standard": "NONE",
                    "mode": "N/A",
                    "reactant_signal": 0,
                    "is_signal": 0,
                    "response_factor": 0,
                    "is_concentration_uM": 0,
                    "concentration_reactant_uM": float('nan'),
                    "products_detected": "N/A",
                    "product_count": 0,
                    "status": "Missing mzML file"
                })
                self.logger.warning(f"Sample {sample_id}: mzML file not found")
                continue
                
            # Find internal standards in this sample - second potential failure point
            present_standards = {}
            for j in range(1, 10):  # Check up to 10 internal standards
                is_col = f"internal-std-{j}"
                if is_col in sample_data.index and not pd.isna(sample_data[is_col]):
                    is_name = sample_data[is_col]
                    present_standards[is_name] = j
        
            if not present_standards:
                # Add entry for sample without internal standards
                conversion_results.append({
                    "sample_id": sample_id,
                    "reactant_name": "UNKNOWN",
                    "reactant_idx": 0,
                    "internal_standard": "NONE",
                    "mode": "N/A",
                    "reactant_signal": 0,
                    "is_signal": 0,
                    "response_factor": 0,
                    "is_concentration_uM": 0,
                    "concentration_reactant_uM": float('nan'),
                    "products_detected": "N/A",
                    "product_count": 0,
                    "status": "No internal standards found"
                })
                self.logger.warning(f"Sample {sample_id}: No internal standards found")
                continue
            
            # Find reactants in this sample - third potential failure point
            reactant_count = self.get_compound_count(sample_data, "reactant")
            
            if reactant_count == 0:
                # Add entry for sample without reactants
                for is_name in present_standards:
                    conversion_results.append({
                        "sample_id": sample_id,
                        "reactant_name": "NONE",
                        "reactant_idx": 0,
                        "internal_standard": is_name,
                        "mode": "N/A",
                        "reactant_signal": 0,
                        "is_signal": 0,
                        "response_factor": 0,
                        "is_concentration_uM": 0,
                        "concentration_reactant_uM": float('nan'),
                        "products_detected": "N/A",
                        "product_count": 0,
                        "status": "No reactants found"
                    })
                self.logger.warning(f"Sample {sample_id}: No reactants found")
                continue
            
            # Get product detection information
            products_found = []
            product_count = self.get_compound_count(sample_data, "product")
            
            # Find all detected products across all modes
            for p_idx in range(1, product_count + 1):
                # Check Positive mode
                p_area_pos = f"product-{p_idx}-EIC-area-Positive"
                if (p_area_pos in sample_data.index and 
                    not pd.isna(sample_data[p_area_pos]) and 
                    float(sample_data[p_area_pos]) > 0):
                    products_found.append(f"Product-{p_idx}-Positive")
                    
                # Check Negative mode
                p_area_neg = f"product-{p_idx}-EIC-area-Negative"
                if (p_area_neg in sample_data.index and 
                    not pd.isna(sample_data[p_area_neg]) and 
                    float(sample_data[p_area_neg]) > 0):
                    products_found.append(f"Product-{p_idx}-Negative")
        
            # Format product information
            products_found_str = ", ".join(products_found) if products_found else "No products detected"
            
            # Process each standard/reactant combination
            for is_name, is_idx in present_standards.items():
                # Lookup this IS in the configuration
                is_config = None
                for is_info in self.config["conversion"]["internal_standards"].get("standards", []):
                    if is_info["name"] == is_name:
                        is_config = is_info
                        break
            
                if is_config is None:
                    # Add entry for IS not in configuration
                    for i in range(1, reactant_count + 1):
                        reactant_col = f"reactant-{i}"
                        reactant_name = sample_data.get(reactant_col, f"Reactant-{i}")
                        conversion_results.append({
                            "sample_id": sample_id,
                            "reactant_name": reactant_name,
                            "reactant_idx": i,
                            "internal_standard": is_name,
                            "mode": "N/A",
                            "reactant_signal": 0,
                            "is_signal": 0,
                            "response_factor": 0,
                            "is_concentration_uM": 0,
                            "concentration_reactant_uM": float('nan'),
                            "products_detected": products_found_str,
                            "product_count": len(products_found),
                            "status": "Internal standard not in configuration"
                        })
                    self.logger.warning(f"Sample {sample_id}: IS {is_name} not found in configuration")
                    continue
                
                # Get IS details
                is_smiles = is_config["smiles"]
                is_modes = is_config.get("modes", ["Positive", "Negative"])
                is_conc = is_config.get("concentration_uM", 
                    self.config["conversion"]["internal_standards"].get("default_concentration_uM", 100.0))
                is_custom_mw = is_config.get("molecular_weight", None)
                
                # Process each reactant
                for i in range(1, reactant_count + 1):
                    reactant_col = f"reactant-{i}"
                    reactant_name = sample_data.get(reactant_col, f"Reactant-{i}")
                    reactant_smiles = sample_data.get(reactant_col)
                    
                    # Check for valid SMILES
                    if not isinstance(reactant_smiles, str) or not reactant_smiles.strip():
                        conversion_results.append({
                            "sample_id": sample_id,
                            "reactant_name": reactant_name,
                            "reactant_idx": i,
                            "internal_standard": is_name,
                            "mode": "N/A",
                            "reactant_signal": 0,
                            "is_signal": 0,
                            "response_factor": 0,
                            "is_concentration_uM": is_conc,
                            "concentration_reactant_uM": float('nan'),
                            "products_detected": products_found_str,
                            "product_count": len(products_found),
                            "status": "Invalid or missing SMILES"
                        })
                        self.logger.warning(f"Sample {sample_id}, Reactant {i}: Invalid or missing SMILES")
                        continue
                    
                    # Standardize SMILES for matching
                    reactant_smiles = desalt_smiles(reactant_smiles)
                    reactant_inchi = standardise_compound(reactant_smiles, logger=self.logger)
                    
                    # Create RF keys - use InChI if available
                    if reactant_inchi:
                        rf_key_pos = f"{reactant_inchi}_{is_name}_Positive"
                        rf_key_neg = f"{reactant_inchi}_{is_name}_Negative"
                    else:
                        rf_key_pos = f"{reactant_name}_{is_name}_Positive"
                        rf_key_neg = f"{reactant_name}_{is_name}_Negative"
                    
                    # Check for response factors in each mode
                    rf_match_pos = rf_key_pos in response_factors
                    rf_match_neg = rf_key_neg in response_factors
                    
                    if not (rf_match_pos or rf_match_neg):
                        # Add entry for missing response factor
                        conversion_results.append({
                            "sample_id": sample_id,
                            "reactant_name": reactant_name,
                            "reactant_idx": i,
                            "internal_standard": is_name,
                            "mode": "BOTH",
                            "reactant_signal": 0,
                            "is_signal": 0,
                            "response_factor": 0,
                            "is_concentration_uM": is_conc,
                            "concentration_reactant_uM": float('nan'),
                            "products_detected": products_found_str,
                            "product_count": len(products_found),
                            "status": "No matching response factor"
                        })
                        self.logger.warning(f"Sample {sample_id}, Reactant {reactant_name}: No matching response factor with {is_name}")
                        continue
                    
                    # Try each mode with a matching RF
                    for mode, rf_key, rf_match in [
                        ("Positive", rf_key_pos, rf_match_pos),
                        ("Negative", rf_key_neg, rf_match_neg)
                    ]:
                        if not rf_match or mode not in is_modes:
                            continue
                        
                        # Extract IS signal
                        is_signal = self._extract_compound_signal(
                            mzML_filepath, is_smiles, mode, custom_mw=is_custom_mw
                        )
                        
                        if is_signal <= 0:
                            # Add entry for no IS signal
                            conversion_results.append({
                                "sample_id": sample_id,
                                "reactant_name": reactant_name,
                                "reactant_idx": i,
                                "internal_standard": is_name,
                                "mode": mode,
                                "reactant_signal": 0,
                                "is_signal": 0,
                                "response_factor": response_factors[rf_key]["response_factor"],
                                "is_concentration_uM": is_conc,
                                "concentration_reactant_uM": float('nan'),
                                "products_detected": products_found_str,
                                "product_count": len(products_found),
                                "status": "No internal standard signal detected"
                            })
                            self.logger.warning(f"Sample {sample_id}, Mode {mode}: No signal for IS {is_name}")
                            continue
                        
                        # Get reactant signal
                        signal_column = f"{reactant_col}-EIC-area-{mode}"
                        
                        if signal_column not in sample_data.index or pd.isna(sample_data[signal_column]):
                            # Add entry for missing reactant signal
                            conversion_results.append({
                                "sample_id": sample_id,
                                "reactant_name": reactant_name,
                                "reactant_idx": i,
                                "internal_standard": is_name,
                                "mode": mode,
                                "reactant_signal": 0,
                                "is_signal": is_signal,
                                "response_factor": response_factors[rf_key]["response_factor"],
                                "is_concentration_uM": is_conc,
                                "concentration_reactant_uM": float('nan'),
                                "products_detected": products_found_str,
                                "product_count": len(products_found),
                                "status": "No reactant signal detected"
                            })
                            self.logger.warning(f"Sample {sample_id}, Mode {mode}: No signal for reactant {reactant_name}")
                            continue
                        
                        # We have all needed data - calculate conversion
                        reactant_signal = float(sample_data[signal_column])
                        response_factor = response_factors[rf_key]["response_factor"]
                        
                        # Calculate concentration
                        concentration_reactant = ((reactant_signal / is_signal) * is_conc) / response_factor
                        
                        # Add successful result
                        conversion_results.append({
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
                            "products_detected": products_found_str,
                            "product_count": len(products_found),
                            "status": "Success"
                        })
                        
                        # Also update batch data
                        col_name = f"reactant-{i}-concentration-{mode}-{is_name}"
                        self.batch_data.loc[idx, col_name] = concentration_reactant
                        
                        self.logger.info(f"Sample {sample_id}: Successfully calculated concentration for {reactant_name} ({concentration_reactant:.2f} µM)")
    
        # Create DataFrame and add status summary
        if conversion_results:
            df = pd.DataFrame(conversion_results)
            
            # Summarize results by status
            status_counts = df.groupby('status').size().to_dict()
            self.logger.info(f"Conversion results by status: {status_counts}")
            self.logger.info(f"Created conversion results DataFrame with {len(df)} rows")
            
            # Calculate success rate
            if 'Success' in status_counts:
                success_rate = status_counts['Success'] / len(df) * 100
                self.logger.info(f"Success rate: {success_rate:.1f}%")
            
            return df
        else:
            self.logger.warning("No conversion results were calculated")
            return pd.DataFrame()
