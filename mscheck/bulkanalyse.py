# %%
from datetime import datetime
import os
import logging
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem

# Local imports
from analyse import AnalyseSpectrum
from report import MSReport
from heatmap import MSHeatmapGenerator
from logging_config import setup_logger, get_logger

# Setup the ROOT logger with a file for this run
log_file = os.path.join("logs", f"mscheck_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
logger = setup_logger("mscheck", level="INFO", log_file=log_file, use_colors=True)

logger.info("Starting MSCheck analysis")

# Set up logging
logger = get_logger(__name__)

class BulkAnalyser:
    """Handles bulk analysis of mass spectrometry data using CSV input"""
    
    def __init__(self, csv_input_path: str, data_dir: str = None, report_dir: str = None):
        """
        Initialize bulk analyzer with paths for input, data, and reports
        
        Args:
            csv_input_path: Path to the CSV file containing sample information
            data_dir: Directory containing mzML files (if None, uses directory of CSV)
            report_dir: Directory to save report files (if None, creates reports subfolder)
        """
        self.csv_input_path = csv_input_path
        self.data_dir = data_dir or os.path.dirname(csv_input_path)
        self.report_dir = report_dir or os.path.join(os.path.dirname(csv_input_path), "reports")
        self.batch_data = None
        self.processed_samples = 0
        self.total_samples = 0
        self.errors = []
        self.logger = logging.getLogger("BulkAnalyser")
        
        self.logger.info(f"Initialized BulkAnalyser with CSV: {csv_input_path}")
        self.logger.info(f"Data directory: {self.data_dir}")
        self.logger.info(f"Report directory: {self.report_dir}")
        
    def load_data(self) -> pd.DataFrame:
        """
        Load data from CSV file
        
        Returns:
            DataFrame with loaded data
        """
        try:
            self.batch_data = pd.read_csv(self.csv_input_path)
            self.total_samples = len(self.batch_data)
            self.logger.info(f"Successfully loaded {self.total_samples} samples from {self.csv_input_path}")
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
        tolerance: int = 1
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
                self.logger.info(f"Processing sample {sample_id} ({batch_index + 1}/{len(self.batch_data)})")
                
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
                        # Get number of compounds of this type
                        no_type_column = f"no-{analysis_type}s"
                        if no_type_column not in batch_row:
                            self.logger.warning(f"'{no_type_column}' column missing for sample {sample_id}")
                            continue
                        
                        no_compounds = int(batch_row[no_type_column])
                        if no_compounds <= 0:
                            continue
                            
                        # Get ion information
                        ions_add_column = f"{analysis_type}-ions-to-add"
                        ions_sub_column = f"{analysis_type}-ions-to-sub"
                        
                        # Get ions with safe parsing
                        analysis_type_ions_to_add = self.parse_ions(batch_row.get(ions_add_column, ""))
                        analysis_type_ions_to_sub = self.parse_ions(batch_row.get(ions_sub_column, ""))
                        
                        # Get tolerance - use default if not specified
                        tolerance_column = f"{analysis_type}-match-tolerance"
                        if tolerance_column in batch_row:
                            analysis_type_match_tolerance = int(batch_row[tolerance_column])
                        else:
                            analysis_type_match_tolerance = tolerance
                            
                        # Process each compound of this type
                        for compound_idx in range(no_compounds):
                            # Get compound SMILES
                            compound_column = f"{analysis_type}-{compound_idx + 1}"
                            if compound_column not in batch_row:
                                self.logger.warning(f"'{compound_column}' column missing for sample {sample_id}")
                                continue
                                
                            analysis_type_smiles = batch_row[compound_column]
                            if not isinstance(analysis_type_smiles, str) or not analysis_type_smiles.strip():
                                self.logger.warning(f"Empty SMILES in {compound_column} for sample {sample_id}")
                                continue
                                
                            # Process each ionization mode
                            for mode in modes:
                                try:
                                    total_analyses_attempted += 1
                                    self.logger.info(f"Analyzing {analysis_type}-{compound_idx + 1} in {mode} mode")
                                    
                                    # Perform analysis
                                    analysis_obj = AnalyseSpectrum(mzMLfilepath=mzML_filepath, mode=mode)
                                    
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
                                        folder=report_subdir,
                                        compound_name=report_name
                                    )
                                    
                                    # Update results in DataFrame using proper column names
                                    signal_column = f"{analysis_type}-{compound_idx + 1}-max-EIC-signal-{mode}"
                                    mz_match_column = f"{analysis_type}-{compound_idx + 1}-max-mz-match-{mode}"
                                    ions_column = f"{analysis_type}-{compound_idx + 1}-ions-matched-{mode}"
                                    eic_area_column = f"{analysis_type}-{compound_idx + 1}-EIC-area-{mode}"
                                    
                                    # Use loc to properly update DataFrame
                                    self.batch_data.loc[batch_index, signal_column] = analysis_obj.analysedata["max_EIC_signal"]
                                    self.batch_data.loc[batch_index, mz_match_column] = str(analysis_obj.analysedata["max_mz_match"])
                                    self.batch_data.loc[batch_index, ions_column] = str(analysis_obj.analysedata["ions"])
                                    self.batch_data.loc[batch_index, eic_area_column] = eic_area
                                    
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
        self.logger.info(f"Analysis complete: {successful_analyses}/{total_analyses_attempted} analyses successful")
        self.logger.info(f"Processed {self.processed_samples}/{self.total_samples} samples")
        
        if self.errors:
            self.logger.warning(f"Encountered {len(self.errors)} errors during processing")
            
        return self.batch_data
        
    def find_mzml_file(self, filename: str) -> Optional[str]:
        """
        Helper method to locate mzML files with flexible path handling
        
        Args:
            filename: Base filename with or without extension
            
        Returns:
            Full path to mzML file if found, None otherwise
        """
        # Try different possible file paths
        possible_paths = [
            os.path.join(self.data_dir, f"{filename}.mzML"),
            os.path.join(self.data_dir, filename),
        ]
        
        # Try with extension if not already present
        if not filename.endswith('.mzML'):
            possible_paths.append(os.path.join(self.data_dir, f"{filename}.mzML"))
        else:
            # If filename already has extension, try without it
            base_name = filename[:-5]
            possible_paths.append(os.path.join(self.data_dir, base_name))
        
        # Check if any path exists
        for path in possible_paths:
            if os.path.exists(path):
                return path
                
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
            return [s.strip() for s in ion_str.split(',') if s.strip()]
        else:
            # Need to split and add brackets if needed
            ions = []
            for ion in ion_str.split(','):
                ion = ion.strip()
                if ion:
                    # Add brackets if needed
                    if not (ion.startswith('[') and ion.endswith(']')):
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
            if os.path.isdir(output_path) or not output_path.lower().endswith('.csv'):
                os.makedirs(output_path, exist_ok=True)
                timestamp = pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')
                output_file = os.path.join(output_path, f"analysis_results_{timestamp}.csv")
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
            "error_messages": self.errors
        }
    
    def visualize_results(
        self, 
        metrics: list = None,
        output_dir: str = None,
        show_plots: bool = False
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
        tolerance: int = 1
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
            analysis_types=analysis_types,
            modes=modes,
            tolerance=tolerance
        )
        
        # Save the results
        results_path = os.path.join(os.path.dirname(self.csv_input_path), "analysis_results.csv")
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
        heatmap_generator = MSHeatmapGenerator(data=self.batch_data, output_dir=self.report_dir)
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
                    "all_mz_data": None
                }
            
            # Get sample chromatogram data if not already loaded
            self._load_sample_chromatogram(samples, sample_id, row)
            
            # Process each compound type
            for compound_type in ["reactant", "product", "internal-std", "intermediate"]:
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
                if 'mz_data' in temp_analyzer.MSdata:
                    samples[sample_id]["all_mz_data"] = temp_analyzer.MSdata["mz_data"]
                    self.logger.info(f"Added complete MZ data to sample {sample_id}")
        except Exception as e:
            self.logger.error(f"Error loading RT/TIC for sample {sample_id}: {str(e)}")
    
    def _extract_compound_data(self, samples, sample_id, row, compound_type):
        """Extract compound data for a specific compound type"""
        # Get number of compounds
        no_type_column = f"no-{compound_type}s"
        if no_type_column not in row:
            return
            
        try:
            no_compounds = int(row[no_type_column])
        except (ValueError, TypeError):
            return
            
        # Process each compound
        for compound_idx in range(no_compounds):
            compound_data = self._extract_single_compound(
                sample_id, row, compound_type, compound_idx)
            
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
                self.logger.warning(f"Could not create molecule from SMILES: {smiles_value}")
                return None
        except Exception as e:
            self.logger.error(f"Error creating molecule: {str(e)}")
            return None
        
        # Get the compound's mass spectrum and other data
        rt_max, intensity_max, eic_data, mz_strongest, all_mz_data = self._get_compound_ms_data(
            sample_id, row, compound_type, compound_idx, smiles_value, mode)
        
        # Create compound entry
        return {
            "compound_name": f"{compound_type.capitalize()} {compound_idx + 1}",
            "compound_type": compound_type,
            "msmode": mode,
            "mol": mol,
            "smiles": smiles_value,
            "rt_max": rt_max,
            "intensity_max": intensity_max,
            "eic_data": eic_data,
            "signal": row.get(signal_column),
            "eic_area": row.get(eic_area_column),
            "max_mz_match": row.get(mz_match_column, "False") == "True",
            "mz_strongest": mz_strongest,
            "all_mz_data": all_mz_data
        }
    
    def _get_compound_ms_data(self, sample_id, row, compound_type, compound_idx, 
                              smiles_value, mode):
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
                tolerance=tolerance
            )
            
            # Extract RT max and intensity max
            if "RT" in temp_analyzer.analysedata and len(temp_analyzer.analysedata["RT"]) > 0:
                rt_values = temp_analyzer.analysedata["RT"][0]
                tic_values = temp_analyzer.analysedata["TIC"][0]
                if len(rt_values) > 0 and len(tic_values) > 0:
                    max_idx = np.argmax(tic_values)
                    rt_max = rt_values[max_idx]
                    intensity_max = tic_values[max_idx]
                    
            # Extract EIC data
            if "EIC_data" in temp_analyzer.analysedata:
                eic_data = temp_analyzer.analysedata["EIC_data"][0]
            
            # Extract mass spectrum data
            if ("mz_strongest" in temp_analyzer.analysedata and 
                len(temp_analyzer.analysedata["mz_strongest"]) > 0):
                mz_masses, mz_intensities, _ = temp_analyzer.analysedata["mz_strongest"][0]
                mz_strongest = (mz_masses, mz_intensities)

            # Get full MZ data
            if hasattr(temp_analyzer, 'MSdata') and 'mz_data' in temp_analyzer.MSdata:
                all_mz_data = temp_analyzer.MSdata['mz_data']
                
        except Exception as e:
            self.logger.error(f"Error extracting MS data: {str(e)}")
            
        return rt_max, intensity_max, eic_data, mz_strongest, all_mz_data
    
    def generate_compound_reports(self, samples=None):
        """Generate multi-compound reports for samples"""
        if samples is None:
            samples = self.extract_sample_data()
            
        report_generator = MSReport(output_dir=os.path.join(self.report_dir, "compound_reports"))
        report_paths = {}
        
        # First pass - collect information about all reports we'll generate
        available_reports = []
        for sample_id, sample_data in samples.items():
            if len(sample_data.get("compounds", [])) > 0 and sample_data.get("RT_values") is not None:
                report_filename = f"Sample_{sample_id}_Analysis.html"
                report_path = os.path.join(self.report_dir, "compound_reports", f"{sample_id}_analysis.html") 
                report_title = f"Sample {sample_id} Analysis"
                
                available_reports.append({
                    "title": report_title,
                    "path": report_filename,
                    "sample_id": sample_id
                })
        
        # Second pass - actually generate the reports with navigation
        for i, report_info in enumerate(available_reports):
            sample_id = report_info["sample_id"]
            sample_data = samples[sample_id]
            
            # Filter valid compounds
            valid_compounds = [c for c in sample_data["compounds"] 
                              if c.get("rt_max") is not None and c.get("mol") is not None]
            
            if not valid_compounds:
                self.logger.warning(f"No compounds with valid RT data for sample {sample_id}")
                continue
                
            self.logger.info(f"Creating report with {len(valid_compounds)} compounds for {sample_id}")
            
            # Generate report
            report_path = report_generator.create_annotated_tic_report(
                RT_values=sample_data["RT_values"],
                TIC_values=sample_data["TIC_values"],
                compounds=valid_compounds,
                report_title=report_info["title"],
                html_output=True,
                svg_layout="triple",
                available_reports=available_reports,
                current_report_index=i
            )
            
            if report_path:
                report_paths[sample_id] = report_path
                self.logger.info(f"Generated report: {report_path}")
            else:
                self.logger.error(f"Failed to generate report for sample {sample_id}")
                
        return report_paths
    
    def run_complete_workflow(
        self,
        analysis_types: List[str] = ["product", "reactant", "internal-std"],
        modes: List[str] = ["Positive"],
        tolerance: int = 1
    ):
        """
        Run the complete analysis workflow
        
        Args:
            analysis_types: List of compound types to analyze (e.g., product, reactant)
            modes: List of ionization modes (e.g., Positive, Negative)
            tolerance: Default mass tolerance
            
        Returns:
            Dictionary mapping sample IDs to their report paths
        """
        self.logger.info("Starting complete analysis workflow...")
        self.logger.info(f"Analysis types: {', '.join(analysis_types)}")
        self.logger.info(f"Modes: {', '.join(modes)}")
        self.logger.info(f"Default tolerance: {tolerance}")
        
        # 1. Process all samples with the specified parameters
        self.process_samples(
            analysis_types=analysis_types,
            modes=modes,
            tolerance=tolerance
        )
        
        # 2. Generate visualizations
        self.generate_visualizations()
        
        # 3. Extract sample data
        samples = self.extract_sample_data()
        
        # 4. Generate compound reports
        report_paths = self.generate_compound_reports(samples)
        
        self.logger.info("\nAnalysis and visualization complete!")
        self.logger.info(f"Reports saved to: {self.report_dir}")
        
        return report_paths


# Set the paths - replace these with your actual paths "
csv_path = "/Users/bvh64415/myrepos/mscheck/tests/testdata/bulk-test/bulk-test-copy.csv"
data_dir = "/Users/bvh64415/myrepos/mscheck/tests/testdata/bulk-test/datafiles/"
report_dir = "/Users/bvh64415/myrepos/mscheck/tests/testdata/bulk-test/reports/"

# Example usage with parameters
logging.info("Starting bulk analysis workflow...")

# Initialize the analyzer
analyzer = BulkAnalyser(
    csv_input_path=csv_path,
    data_dir=data_dir,
    report_dir=report_dir
)

# Run the complete workflow with custom parameters
report_paths = analyzer.run_complete_workflow(
    analysis_types=["product", "reactant", "internal-std", "intermediate"],
    modes=["Positive", "Negative"],
    tolerance=1
)

logging.info("\nAnalysis and visualization complete!")
logging.info(f"Reports saved to: {report_dir}")
# %%