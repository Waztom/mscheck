# %%
import os
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Union
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import LogNorm

from analyse import AnalyseSpectrum

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
        
    def load_data(self) -> pd.DataFrame:
        """
        Load data from CSV file
        
        Returns:
            DataFrame with loaded data
        """
        try:
            self.batch_data = pd.read_csv(self.csv_input_path)
            self.total_samples = len(self.batch_data)
            print(f"Successfully loaded {self.total_samples} samples from {self.csv_input_path}")
            return self.batch_data
        except Exception as e:
            error_msg = f"Error loading CSV file: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)
            self.batch_data = pd.DataFrame()
            return self.batch_data
            
    def analyse_batch(
        self, 
        analysis_types: List[str] = ["product", "intermediate", "reactant"], 
        modes: List[str] = ["Positive"], 
        default_tolerance: int = 1
    ) -> pd.DataFrame:
        """
        Process all samples in the batch
        
        Args:
            analysis_types: List of analysis types to perform (e.g., product, reactant)
            modes: List of ionization modes to analyze
            default_tolerance: Default mass tolerance if not specified in CSV
            
        Returns:
            DataFrame with analysis results
        """
        if self.batch_data is None:
            self.load_data()
            
        if len(self.batch_data) == 0:
            print("No data to analyze")
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
                print(f"Processing sample {sample_id} ({batch_index + 1}/{len(self.batch_data)})")
                
                # Find mzML file
                mzML_filepath = self.find_mzml_file(mzML_filename)
                if mzML_filepath is None:
                    error_msg = f"Warning: mzML file not found for sample {sample_id}: {mzML_filename}"
                    self.errors.append(error_msg)
                    print(error_msg)
                    continue
                
                # Process each analysis type
                for analysis_type in analysis_types:
                    try:
                        # Get number of compounds of this type
                        no_type_column = f"no-{analysis_type}s"
                        if no_type_column not in batch_row:
                            print(f"Warning: '{no_type_column}' column missing for sample {sample_id}")
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
                            analysis_type_match_tolerance = default_tolerance
                            
                        # Process each compound of this type
                        for compound_idx in range(no_compounds):
                            # Get compound SMILES
                            compound_column = f"{analysis_type}-{compound_idx + 1}"
                            if compound_column not in batch_row:
                                print(f"Warning: '{compound_column}' column missing for sample {sample_id}")
                                continue
                                
                            analysis_type_smiles = batch_row[compound_column]
                            if not isinstance(analysis_type_smiles, str) or not analysis_type_smiles.strip():
                                print(f"Warning: Empty SMILES in {compound_column} for sample {sample_id}")
                                continue
                                
                            # Process each ionization mode
                            for mode in modes:
                                try:
                                    total_analyses_attempted += 1
                                    print(f"  Analyzing {analysis_type}-{compound_idx + 1} in {mode} mode")
                                    
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
                                    print(f"    EIC Area: {eic_area:.2f}")
                                    
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
                                    print(error_msg)
                    except Exception as e:
                        error_msg = f"Error processing {analysis_type} for sample {sample_id}: {str(e)}"
                        self.errors.append(error_msg)
                        print(error_msg)
                
                self.processed_samples += 1
                
            except Exception as e:
                error_msg = f"Error processing sample {batch_index}: {str(e)}"
                self.errors.append(error_msg)
                print(error_msg)
        
        # Print summary
        print(f"Analysis complete: {successful_analyses}/{total_analyses_attempted} analyses successful")
        print(f"Processed {self.processed_samples}/{self.total_samples} samples")
        
        if self.errors:
            print(f"Encountered {len(self.errors)} errors during processing")
            
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
            print("No data to save")
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
            print(f"Results saved to {output_file}")
            return True
            
        except Exception as e:
            error_msg = f"Error saving results: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)
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
    
    def generate_heatmap(
        self,
        signal_column: str,
        output_dir: str = None,
        plate_column: str = "plate-ID",
        well_column: str = "well-ID",
        title: str = None,
        colormap: str = "viridis",
        save_format: str = "png",
        dpi: int = 300,
        show_plot: bool = True
    ) -> None:
        """
        Generate a heatmap visualization of signals across plates and wells
        
        Args:
            signal_column: Column name containing the signal to visualize
            output_dir: Directory to save heatmap figures (default: report_dir/heatmaps)
            plate_column: Column containing plate IDs
            well_column: Column containing well IDs (format: A1, B2, etc.)
            title: Custom title for the plot (default: derived from signal_column)
            colormap: Matplotlib colormap to use
            save_format: Format for saving figures ('png', 'pdf', 'svg')
            dpi: Resolution for saving figures
            show_plot: Whether to display the plot (set to False for batch processing)
        """
        if self.batch_data is None or len(self.batch_data) == 0:
            print("No data to visualize")
            return
        
        if output_dir is None:
            output_dir = os.path.join(self.report_dir, "heatmaps")
        
        os.makedirs(output_dir, exist_ok=True)
        
        if title is None:
            title = f"Heatmap of {signal_column}"
        
        # Check if required columns exist for plate layout
        has_plate_layout = plate_column in self.batch_data.columns and well_column in self.batch_data.columns
        
        if has_plate_layout:
            # Generate plate-based heatmap
            self._generate_plate_heatmap(
                signal_column=signal_column,
                plate_column=plate_column,
                well_column=well_column,
                output_dir=output_dir,
                title=title,
                colormap=colormap,
                save_format=save_format,
                dpi=dpi,
                show_plot=show_plot
            )
        else:
            # Generate simple heatmap with just the values
            self._generate_simple_heatmap(
                signal_column=signal_column,
                output_dir=output_dir,
                title=title,
                colormap=colormap,
                save_format=save_format,
                dpi=dpi,
                show_plot=show_plot
            )

    def _generate_simple_heatmap(
        self,
        signal_column: str,
        output_dir: str,
        title: str,
        colormap: str,
        save_format: str,
        dpi: int,
        show_plot: bool
    ) -> None:
        """Generate a simple heatmap without plate/well information"""
        try:
            if signal_column not in self.batch_data.columns:
                print(f"Warning: Column '{signal_column}' not found in data")
                return
                
            # Extract numeric data and remove NaNs
            data_values = pd.to_numeric(self.batch_data[signal_column], errors='coerce')
            valid_data = data_values.dropna()
            
            if len(valid_data) == 0:
                print(f"Warning: No valid numeric data in column '{signal_column}'")
                return
            
            # Create a matrix suitable for heatmap display
            # We'll make it as square as possible
            size = int(np.ceil(np.sqrt(len(valid_data))))
            data_matrix = np.zeros((size, size))
            data_matrix.fill(np.nan)
            
            # Fill the matrix with values
            for i, val in enumerate(valid_data):
                row = i // size
                col = i % size
                if row < size and col < size:
                    data_matrix[row, col] = val
            
            # Create figure with constrained_layout
            fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
            
            # Determine if log scale would be appropriate
            non_nan_values = data_matrix[~np.isnan(data_matrix)]
            if len(non_nan_values) > 0:
                min_val = np.nanmin(non_nan_values)
                max_val = np.nanmax(non_nan_values)
                
                if min_val > 0 and max_val / min_val > 100:
                    # Use log scale for large dynamic range
                    norm = LogNorm(vmin=max(min_val, 0.1), vmax=max_val)
                    cmap_label = f"{signal_column} (log scale)"
                else:
                    # Use linear scale
                    norm = None
                    cmap_label = signal_column
            else:
                norm = None
                cmap_label = signal_column
            
            # Create the heatmap with explicit axes
            sns.heatmap(
                data_matrix, 
                cmap=colormap,
                annot=True,
                fmt=".2g",
                linewidths=0.5,
                cbar_kws={'label': cmap_label},
                norm=norm,
                mask=np.isnan(data_matrix),
                ax=ax  # Explicitly pass the axes
            )
            
            ax.set_title(title)
            ax.set_xlabel("Column Index")
            ax.set_ylabel("Row Index")
            
            # Create a sanitized filename
            safe_signal_column = signal_column.replace('-', '_').replace(' ', '_')
            
            # Save figure
            output_file = os.path.join(output_dir, f"heatmap_{safe_signal_column}.{save_format}")
            fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
            print(f"Generated simple heatmap: {output_file}")
            
            if show_plot:
                plt.show()
            else:
                plt.close(fig)
                
        except Exception as e:
            error_msg = f"Error generating simple heatmap: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)

    def _generate_plate_heatmap(
        self,
        signal_column: str,
        plate_column: str,
        well_column: str,
        output_dir: str,
        title: str,
        colormap: str,
        save_format: str,
        dpi: int,
        show_plot: bool
    ) -> None:
        """Generate plate-organized heatmaps"""
        try:
            # Check if required columns exist
            required_cols = [plate_column, well_column, signal_column]
            missing_cols = [col for col in required_cols if col not in self.batch_data.columns]
            if missing_cols:
                print(f"Warning: Missing required columns: {', '.join(missing_cols)}")
                return
            
            # Get unique plates
            unique_plates = self.batch_data[plate_column].unique()
            
            for plate_id in unique_plates:
                # Filter data for this plate
                plate_data = self.batch_data[self.batch_data[plate_column] == plate_id].copy()
                
                if len(plate_data) == 0:
                    continue
                
                # Determine plate format based on well IDs
                # Get the highest row and column to determine plate format
                highest_well = max(plate_data[well_column].dropna(), key=lambda x: str(x) if isinstance(x, str) else "")
                if isinstance(highest_well, str) and len(highest_well) >= 2:
                    highest_row_letter = highest_well[0].upper()
                    highest_row = ord(highest_row_letter) - ord('A')
                    
                    # Extract the numeric part - could be more than one digit (e.g., A12, P24)
                    col_str = ''.join(c for c in highest_well[1:] if c.isdigit())
                    highest_col = int(col_str) - 1 if col_str else 0
                    
                    # Determine plate format
                    if highest_row >= 15 or highest_col >= 23:  # 384-well plate (P24 is the highest well)
                        n_rows = 16
                        n_cols = 24
                        plate_type = "384-well"
                    else:  # 96-well plate (H12 is the highest well)
                        n_rows = 8
                        n_cols = 12
                        plate_type = "96-well"
                else:
                    # Default to 96-well if can't determine
                    n_rows = 8
                    n_cols = 12
                    plate_type = "96-well (default)"
                
                # Convert well IDs to row/column indices (A1 -> 0,0; B2 -> 1,1; etc.)
                def parse_well(well):
                    if not isinstance(well, str) or len(well) < 2:
                        return -1, -1
                    
                    row_idx = ord(well[0].upper()) - ord('A')
                    
                    # Extract digits from well ID
                    col_str = ''.join(c for c in well[1:] if c.isdigit())
                    col_idx = int(col_str) - 1 if col_str else -1
                    
                    return row_idx, col_idx
                
                # Create empty heatmap matrix with proper dimensions
                heatmap_data = np.zeros((n_rows, n_cols))
                heatmap_data[:] = np.nan  # Set all to NaN initially
                
                # Fill in the data
                for _, row in plate_data.iterrows():
                    if pd.notna(row[signal_column]) and pd.notna(row[well_column]):
                        try:
                            well_id = str(row[well_column])
                            r_idx, c_idx = parse_well(well_id)
                            
                            if r_idx >= 0 and c_idx >= 0 and r_idx < n_rows and c_idx < n_cols:
                                value = float(row[signal_column])
                                heatmap_data[r_idx, c_idx] = value
                        except (ValueError, TypeError):
                            continue
                
                # Create figure with proper aspect ratio for plate layout
                fig_width = max(10, n_cols * 0.6)
                fig_height = max(8, n_rows * 0.6)
                
                # Use constrained_layout instead of tight_layout to avoid the warning
                fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
                
                # Determine if log scale would be appropriate
                non_nan_values = heatmap_data[~np.isnan(heatmap_data)]
                if len(non_nan_values) > 0:
                    min_val = np.nanmin(non_nan_values)
                    max_val = np.nanmax(non_nan_values)
                    
                    if min_val > 0 and max_val / min_val > 100:
                        # Use log scale for large dynamic range
                        norm = LogNorm(vmin=max(min_val, 0.1), vmax=max_val)
                        cmap_label = f"{signal_column} (log scale)"
                    else:
                        # Use linear scale
                        norm = None
                        cmap_label = signal_column
                else:
                    norm = None
                    cmap_label = signal_column
                
                # Determine if we should show values in cells
                # For 384-well plates, the cells might be too small for annotations
                show_values = n_rows <= 8 and n_cols <= 12
                
                # Create the heatmap with explicit axes
                sns.heatmap(
                    heatmap_data, 
                    cmap=colormap,
                    annot=show_values,
                    annot_kws={"size": 8},
                    fmt=".2g",
                    linewidths=0.5,
                    cbar_kws={'label': cmap_label},
                    norm=norm,
                    mask=np.isnan(heatmap_data),
                    ax=ax  # Explicitly pass the axes
                )
                
                # Add well labels
                row_labels = [chr(i + ord('A')) for i in range(n_rows)]
                col_labels = [str(i + 1) for i in range(n_cols)]
                
                # Set tick frequency based on plate size
                row_tick_freq = 1 if n_rows <= 8 else 2
                col_tick_freq = 1 if n_cols <= 12 else 2
                
                # Show only some of the labels for larger plates
                ax.set_yticks(np.arange(n_rows)[::row_tick_freq] + 0.5)
                ax.set_yticklabels(row_labels[::row_tick_freq], rotation=0)
                
                ax.set_xticks(np.arange(n_cols)[::col_tick_freq] + 0.5)
                ax.set_xticklabels(col_labels[::col_tick_freq])
                
                # Set titles and labels
                ax.set_title(f"{title}\nPlate: {plate_id} ({plate_type})")
                
                # Create a sanitized filename
                safe_signal_column = signal_column.replace('-', '_').replace(' ', '_')
                safe_plate_id = str(plate_id).replace('-', '_').replace(' ', '_')
                
                # Save figure
                output_file = os.path.join(output_dir, f"heatmap_{safe_signal_column}_{safe_plate_id}.{save_format}")
                fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
                print(f"Generated heatmap for plate {plate_id}: {output_file}")
                
                if show_plot:
                    plt.show()
                else:
                    plt.close(fig)
        
        except Exception as e:
            error_msg = f"Error generating plate heatmap: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)
            # Print traceback for easier debugging
            import traceback
            traceback.print_exc()

    def _generate_plate_comparison_chart(
        self, 
        signal_column: str, 
        plate_column: str,
        output_dir: str,
        title: str,
        save_format: str,
        dpi: int,
        show_plot: bool
    ) -> None:
        """Generate a comparison chart for multiple plates"""
        try:
            # Get unique plates
            unique_plates = self.batch_data[plate_column].unique()
            
            # Calculate statistics per plate
            plate_stats = []
            for plate in unique_plates:
                plate_data = self.batch_data[self.batch_data[plate_column] == plate]
                values = pd.to_numeric(plate_data[signal_column], errors='coerce').dropna()
                
                if len(values) > 0:
                    plate_stats.append({
                        'plate': plate,
                        'mean': values.mean(),
                        'median': values.median(),
                        'max': values.max(),
                        'min': values.min(),
                        'std': values.std(),
                        'count': len(values)
                    })
            
            if not plate_stats:
                return
                
            # Create a DataFrame for easier plotting
            stats_df = pd.DataFrame(plate_stats)
            
            # Create comparison plot with constrained_layout
            fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
            
            # Plot bars without error bars
            bars = sns.barplot(x='plate', y='mean', data=stats_df, ax=ax)
            
            # Add error bars manually
            for i, row in stats_df.iterrows():
                ax.errorbar(x=i, y=row['mean'], yerr=row['std'], 
                            fmt='none', color='black', capsize=4, capthick=1.5)
            
            # Add value labels on top of bars
            for i, row in stats_df.iterrows():
                ax.text(i, row['mean'] + row['std'] * 0.5, f"{row['mean']:.2f}", 
                        ha='center', va='bottom', fontweight='bold')
            
            ax.set_title(f"Comparison of {title} Across Plates")
            ax.set_ylabel(f"Mean {signal_column}")
            ax.set_xlabel("Plate ID")
            
            # Create a sanitized filename
            safe_signal_column = signal_column.replace('-', '_').replace(' ', '_')
            
            # Save figure
            output_file = os.path.join(output_dir, f"plate_comparison_{safe_signal_column}.{save_format}")
            fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
            print(f"Generated plate comparison chart: {output_file}")
            
            if show_plot:
                plt.show()
            else:
                plt.close(fig)
        
        except Exception as e:
            error_msg = f"Error generating plate comparison: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)
            # Add traceback for debugging
            import traceback
            traceback.print_exc()


# Set the paths - replace these with your actual paths "
csv_path = "/Users/bvh64415/myrepos/mscheck/tests/testdata/bulk-test/bulk-test.csv"
data_dir = "/Users/bvh64415/myrepos/mscheck/tests/testdata/bulk-test/datafiles/"
report_dir = "/Users/bvh64415/myrepos/mscheck/tests/testdata/bulk-test/reports/"

print("Starting bulk analysis workflow...")

# Initialize the analyzer
analyzer = BulkAnalyser(
    csv_input_path=csv_path,
    data_dir=data_dir,
    report_dir=report_dir
)

# Load the CSV data
print("Loading data...")
analyzer.load_data()

# Run the analysis
print("Running batch analysis...")
analyzer.analyse_batch(
    analysis_types=["product", "reactant", "internal-std"],
    modes=["Positive"],
    default_tolerance=1
)

# Save the analysis results
results_path = os.path.join(os.path.dirname(csv_path), "analysis_results.csv")
print(f"Saving results to {results_path}...")
analyzer.save_results(results_path)

# Get analysis summary
summary = analyzer.get_summary()
print("\nAnalysis Summary:")
print(f"Total samples: {summary['total_samples']}")
print(f"Processed samples: {summary['processed_samples']}")
print(f"Errors encountered: {summary['errors']}")

# Generate heatmaps for various metrics
print("\nGenerating heatmaps...")

# Add this before your loop that generates heatmaps
import matplotlib
# Use a different backend if needed
matplotlib.use('Agg')  # Use non-interactive backend for better compatibility

# Update your output path to ensure it's absolute
heatmap_dir = os.path.abspath(os.path.join(report_dir, "platemaps"))
os.makedirs(heatmap_dir, exist_ok=True)
print(f"✓ Created output directory at: {heatmap_dir}")

# Test ability to write to the directory
try:
    test_file = os.path.join(heatmap_dir, "test_write.txt")
    with open(test_file, 'w') as f:
        f.write("Test")
    os.remove(test_file)
    print("✓ Successfully verified write permission to heatmap directory")
except Exception as e:
    print(f"⚠️ Warning: May not be able to write to heatmap directory: {str(e)}")

# When setting up the metrics, use more specific column references:
metrics = [
    # Format: (column_name_prefix, title_prefix, colormap)
    ("product-1-max-EIC-signal-Positive", "Product 1 Signal", "rocket"),
    ("product-1-EIC-area-Positive", "Product 1 Area", "plasma"),
    ("reactant-1-max-EIC-signal-Positive", "Reactant 1 Signal", "inferno"),
    ("reactant-1-EIC-area-Positive", "Reactant 1 Area", "magma")
]

# Before generating heatmaps, print column names for verification:
print("\nAvailable columns in data:")
for col in analyzer.batch_data.columns:
    print(f" - {col}")

# Generate heatmaps for all metrics in all modes
for column_prefix, title_prefix, cmap in metrics:
    column = column_prefix
    
    # Check if column exists in data
    if column in analyzer.batch_data.columns:
        print(f"Generating heatmap for {column}...")
        
        # Generate heatmap
        analyzer.generate_heatmap(
            signal_column=column,
            output_dir=heatmap_dir,
            title=f"{title_prefix}",
            colormap=cmap,
            save_format="png",
            dpi=300,
            show_plot=False  # Set to True to display plots interactively
        )

# Calculate the ratio of product to reactant (conversion)
if "product-1-EIC-area-Positive" in analyzer.batch_data.columns and \
    "reactant-1-EIC-area-Positive" in analyzer.batch_data.columns:
    
    print("\nCalculating conversion ratios...")
    
    # Create conversion ratio column
    analyzer.batch_data["conversion-ratio"] = pd.to_numeric(
        analyzer.batch_data["product-1-EIC-area-Positive"], 
        errors="coerce"
    ) / pd.to_numeric(
        analyzer.batch_data["reactant-1-EIC-area-Positive"], 
        errors="coerce"
    )
    
    # Generate heatmap for conversion ratio
    analyzer.generate_heatmap(
        signal_column="conversion-ratio",
        output_dir=heatmap_dir,
        title="Conversion Ratio (Product/Reactant)",
        colormap="coolwarm",
        save_format="png",
        dpi=300,
        show_plot=False
    )
    
    # Save updated results with conversion data
    analyzer.save_results(os.path.join(os.path.dirname(csv_path), "analysis_with_conversion.csv"))

print("\nGenerating plate comparison charts...")

# Define metrics to compare across plates (same ones used for heatmaps)
comparison_metrics = [
    # Format: (column_name, title, colormap)
    ("product-1-max-EIC-signal-Positive", "Product 1 Signal", "rocket"),
    ("product-1-EIC-area-Positive", "Product 1 Area", "plasma"),
    ("reactant-1-max-EIC-signal-Positive", "Reactant 1 Signal", "inferno"),
    ("reactant-1-EIC-area-Positive", "Reactant 1 Area", "magma"),
    ("conversion-ratio", "Conversion Ratio (Product/Reactant)", "coolwarm")
]

# Create a dedicated directory for comparison charts
comparison_dir = os.path.join(report_dir, "plate_comparisons")
os.makedirs(comparison_dir, exist_ok=True)
print(f"✓ Created directory for plate comparisons: {comparison_dir}")

# Generate comparison charts directly (without needing to create heatmaps first)
for column, title, cmap in comparison_metrics:
    if column in analyzer.batch_data.columns:
        print(f"Generating plate comparison for {title}...")
        
        # Check for plate ID column
        if "plate-ID" in analyzer.batch_data.columns:
            try:
                # Get unique plates and calculate stats
                unique_plates = analyzer.batch_data["plate-ID"].unique()
                plate_stats = []
                
                for plate in unique_plates:
                    plate_data = analyzer.batch_data[analyzer.batch_data["plate-ID"] == plate]
                    values = pd.to_numeric(plate_data[column], errors='coerce').dropna()
                    
                    if len(values) > 0:
                        plate_stats.append({
                            'plate': plate,
                            'mean': values.mean(),
                            'median': values.median(),
                            'max': values.max(),
                            'min': values.min(),
                            'std': values.std(),
                            'count': len(values)
                        })
                
                if plate_stats:
                    # Create a DataFrame for plotting
                    stats_df = pd.DataFrame(plate_stats)
                    
                    # Create comparison plot with constrained_layout
                    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
                    
                    # Plot bars without error bars in the barplot call
                    bars = sns.barplot(x='plate', y='mean', data=stats_df, ax=ax, palette=cmap)
                    
                    # Add error bars manually
                    for i, row in stats_df.iterrows():
                        ax.errorbar(x=i, y=row['mean'], yerr=row['std'], 
                                    fmt='none', color='black', capsize=4, capthick=1.5)
                    
                    # Add value labels on top of bars
                    for i, row in stats_df.iterrows():
                        ax.text(i, row['mean'] + row['std'] * 0.5, f"{row['mean']:.2f}", 
                                ha='center', va='bottom', fontweight='bold')
                    
                    # Add count of samples below each bar
                    for i, row in stats_df.iterrows():
                        ax.text(i, -0.05 * ax.get_ylim()[1], f"n={row['count']}", 
                                ha='center', va='top', fontsize=8)
                    
                    # Add statistical summary as text box
                    summary_text = (
                        f"Overall Mean: {stats_df['mean'].mean():.2f}\n"
                        f"Overall Std: {stats_df['mean'].std():.2f}\n"
                        f"Min: {stats_df['min'].min():.2f}\n"
                        f"Max: {stats_df['max'].max():.2f}\n"
                        f"Total Samples: {stats_df['count'].sum()}"
                    )
                    
                    # Position text box in upper right
                    ax.text(0.95, 0.95, summary_text, transform=ax.transAxes,
                            verticalalignment='top', horizontalalignment='right',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                    
                    ax.set_title(f"Comparison of {title} Across Plates")
                    ax.set_ylabel(column)
                    ax.set_xlabel("Plate ID")
                    
                    # Rotate x-tick labels if needed
                    if len(unique_plates) > 5 or any(len(str(p)) > 10 for p in unique_plates):
                        plt.xticks(rotation=45, ha='right')
                    
                    # Create a sanitized filename
                    safe_column = column.replace('-', '_').replace(' ', '_')
                    
                    # Save the figure
                    output_file = os.path.join(comparison_dir, f"plate_comparison_{safe_column}.png")
                    fig.savefig(output_file, dpi=300, bbox_inches='tight')
                    print(f"✓ Saved plate comparison chart: {output_file}")
                    
                    # Close the figure to free memory
                    plt.close(fig)
                else:
                    print(f"No valid data for {column} across plates")
            except Exception as e:
                print(f"Error creating plate comparison for {column}: {str(e)}")
                import traceback
                traceback.print_exc()
        else:
            print("Cannot create plate comparison - no plate-ID column found")
    else:
        print(f"Column {column} not found in data")

print("\nPlate comparison charts generated!")
print("\nAnalysis and visualization complete!")
print(f"Reports saved to: {report_dir}")
print(f"Heatmaps saved to: {heatmap_dir}")
print(f"Plate comparisons saved to: {comparison_dir}")
# %%