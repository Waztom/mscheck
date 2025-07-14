import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
import traceback
from logging_config import get_logger

# Get logger for this module
logger = get_logger(__name__)


class MSHeatmapGenerator:
    """Handles generation of heatmaps and plate comparison charts for MS analysis data"""

    def __init__(self, data: pd.DataFrame, output_dir: str = None):
        """
        Initialize heatmap generator

        Args:
            data: DataFrame containing analysis results
            output_dir: Base directory for saving visualizations
        """
        self.data = data.copy()
        self.output_dir = output_dir
        self.errors = []
        self.logger = logger

        # Log initialization
        self.logger.info(
            f"Initialized MSHeatmapGenerator with output directory: {output_dir}"
        )

        # Create output directories
        if self.output_dir:
            self.heatmap_dir = os.path.join(self.output_dir, "heatmaps")
            self.comparison_dir = os.path.join(self.output_dir, "comparisons")
            os.makedirs(self.heatmap_dir, exist_ok=True)
            os.makedirs(self.comparison_dir, exist_ok=True)
            self.logger.debug(
                f"Created output directories: {self.heatmap_dir}, {self.comparison_dir}"
            )

    def generate_heatmap(
        self,
        signal_column: str,
        plate_column: str = "plate-ID",
        well_column: str = "well-ID",
        title: str = None,
        colormap: str = "viridis",
        output_dir: str = None,
        save_format: str = "png",
        dpi: int = 300,
        show_plot: bool = True,
    ) -> None:
        """
        Generate a heatmap visualization of signals across plates and wells

        Args:
            signal_column: Column name containing the signal to visualize
            plate_column: Column containing plate IDs
            well_column: Column containing well IDs (format: A1, B2, etc.)
            title: Custom title for the plot (default: derived from signal_column)
            colormap: Matplotlib colormap to use
            output_dir: Custom output directory (overrides default)
            save_format: Format for saving figures ('png', 'pdf', 'svg')
            dpi: Resolution for saving figures
            show_plot: Whether to display the plot
        """
        if self.data is None or len(self.data) == 0:
            print("No data to visualize")
            return

        if output_dir is None:
            output_dir = self.heatmap_dir

        os.makedirs(output_dir, exist_ok=True)

        if title is None:
            title = f"Heatmap of {signal_column}"

        # Check if required columns exist for plate layout
        has_plate_layout = (
            plate_column in self.data.columns and well_column in self.data.columns
        )

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
                show_plot=show_plot,
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
                show_plot=show_plot,
            )

    def _generate_simple_heatmap(
        self,
        signal_column: str,
        output_dir: str,
        title: str,
        colormap: str,
        save_format: str,
        dpi: int,
        show_plot: bool,
    ) -> None:
        """Generate a simple heatmap without plate/well information"""
        try:
            if signal_column not in self.data.columns:
                print(f"Warning: Column '{signal_column}' not found in data")
                return

            # Extract numeric data and remove NaNs
            data_values = pd.to_numeric(self.data[signal_column], errors="coerce")
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
                cbar_kws={"label": cmap_label},
                norm=norm,
                mask=np.isnan(data_matrix),
                ax=ax,  # Explicitly pass the axes
            )

            ax.set_title(title)
            ax.set_xlabel("Column Index")
            ax.set_ylabel("Row Index")

            # Create a sanitized filename
            safe_signal_column = signal_column.replace("-", "_").replace(" ", "_")

            # Save figure
            output_file = os.path.join(
                output_dir, f"heatmap_{safe_signal_column}.{save_format}"
            )
            fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
            print(f"Generated simple heatmap: {output_file}")

            if show_plot:
                plt.show()
            else:
                plt.close(fig)

        except Exception as e:
            error_msg = f"Error generating simple heatmap: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)
            traceback.print_exc()

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
        show_plot: bool,
    ) -> None:
        """Generate plate-organized heatmaps"""
        try:
            # Check if required columns exist
            required_cols = [plate_column, well_column, signal_column]
            missing_cols = [
                col for col in required_cols if col not in self.data.columns
            ]
            if missing_cols:
                print(f"Warning: Missing required columns: {', '.join(missing_cols)}")
                return

            # Get unique plates
            unique_plates = self.data[plate_column].unique()

            for plate_id in unique_plates:
                # Filter data for this plate
                plate_data = self.data[self.data[plate_column] == plate_id].copy()

                if len(plate_data) == 0:
                    continue

                # Determine plate format based on well IDs
                # Get the highest row and column to determine plate format
                highest_well = max(
                    plate_data[well_column].dropna(),
                    key=lambda x: str(x) if isinstance(x, str) else "",
                )
                if isinstance(highest_well, str) and len(highest_well) >= 2:
                    highest_row_letter = highest_well[0].upper()
                    highest_row = ord(highest_row_letter) - ord("A")

                    # Extract the numeric part - could be more than one digit (e.g., A12, P24)
                    col_str = "".join(c for c in highest_well[1:] if c.isdigit())
                    highest_col = int(col_str) - 1 if col_str else 0

                    # Determine plate format
                    if (
                        highest_row >= 15 or highest_col >= 23
                    ):  # 384-well plate (P24 is the highest well)
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

                # Parse well IDs into row/column indices
                def parse_well(well):
                    if not isinstance(well, str) or len(well) < 2:
                        return -1, -1

                    row_idx = ord(well[0].upper()) - ord("A")

                    # Extract digits from well ID
                    col_str = "".join(c for c in well[1:] if c.isdigit())
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

                            if (
                                r_idx >= 0
                                and c_idx >= 0
                                and r_idx < n_rows
                                and c_idx < n_cols
                            ):
                                value = float(row[signal_column])
                                heatmap_data[r_idx, c_idx] = value
                        except (ValueError, TypeError):
                            continue

                # Create figure with proper aspect ratio for plate layout
                fig_width = max(10, n_cols * 0.6)
                fig_height = max(8, n_rows * 0.6)
                fig, ax = plt.subplots(
                    figsize=(fig_width, fig_height), constrained_layout=True
                )

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

                # Create the heatmap
                sns.heatmap(
                    heatmap_data,
                    cmap=colormap,
                    annot=show_values,  # Only show values for 96-well plates
                    annot_kws={"size": 8},  # Smaller font for values
                    fmt=".2g",
                    linewidths=0.5,
                    cbar_kws={"label": cmap_label},
                    norm=norm,
                    mask=np.isnan(heatmap_data),
                    ax=ax,
                )

                # Add well labels - more specific formatting for plate layouts
                row_labels = [chr(i + ord("A")) for i in range(n_rows)]
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
                safe_signal_column = self.sanitize_filename(signal_column)
                safe_plate_id = self.sanitize_filename(str(plate_id))

                # Save figure
                output_file = os.path.join(
                    output_dir,
                    f"heatmap_{safe_signal_column}_{safe_plate_id}.{save_format}",
                )
                fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
                print(f"Generated heatmap for plate {plate_id}: {output_file}")

                if show_plot:
                    plt.show()
                else:
                    plt.close(fig)

        except Exception as e:
            error_msg = f"Error generating plate heatmap: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)
            traceback.print_exc()

    def generate_plate_comparison(
        self,
        signal_column: str,
        plate_column: str = "plate-ID",
        title: str = None,
        colormap: str = "viridis",
        output_dir: str = None,
        save_format: str = "png",
        dpi: int = 300,
        show_plot: bool = True,
    ) -> None:
        """Generate comparison charts across plates"""
        try:
            if output_dir is None:
                output_dir = self.comparison_dir

            os.makedirs(output_dir, exist_ok=True)

            if title is None:
                title = signal_column

            # Check if columns exist
            if plate_column not in self.data.columns:
                print(f"Warning: Plate column '{plate_column}' not found")
                return

            if signal_column not in self.data.columns:
                print(f"Warning: Signal column '{signal_column}' not found")
                return

            # Get unique plates
            unique_plates = self.data[plate_column].unique()

            if len(unique_plates) == 0:
                print("No plates found in data")
                return

            # Calculate statistics per plate
            plate_stats = []
            for plate in unique_plates:
                plate_data = self.data[self.data[plate_column] == plate]
                values = pd.to_numeric(
                    plate_data[signal_column], errors="coerce"
                ).dropna()

                if len(values) > 0:
                    plate_stats.append(
                        {
                            "plate": plate,
                            "mean": values.mean(),
                            "median": values.median(),
                            "max": values.max(),
                            "min": values.min(),
                            "std": values.std(),
                            "count": len(values),
                        }
                    )

            if not plate_stats:
                print(f"No valid data for {signal_column} across plates")
                return

            # Create a DataFrame for easier plotting
            stats_df = pd.DataFrame(plate_stats)

            # Create comparison plot with constrained_layout
            fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

            # Plot bars without error bars
            bars = sns.barplot(
                x="plate", y="mean", data=stats_df, ax=ax, palette=colormap
            )

            # Add error bars manually
            for i, row in stats_df.iterrows():
                ax.errorbar(
                    x=i,
                    y=row["mean"],
                    yerr=row["std"],
                    fmt="none",
                    color="black",
                    capsize=4,
                    capthick=1.5,
                )

            # Add value labels on top of bars
            for i, row in stats_df.iterrows():
                ax.text(
                    i,
                    row["mean"] + row["std"] * 0.5,
                    f"{row['mean']:.2f}",
                    ha="center",
                    va="bottom",
                    fontweight="bold",
                )

            # Add count of samples below each bar
            for i, row in stats_df.iterrows():
                ax.text(
                    i,
                    -0.05 * ax.get_ylim()[1],
                    f"n={row['count']}",
                    ha="center",
                    va="top",
                    fontsize=8,
                )

            # Add statistical summary as text box
            summary_text = (
                f"Overall Mean: {stats_df['mean'].mean():.2f}\n"
                f"Overall Std: {stats_df['mean'].std():.2f}\n"
                f"Min: {stats_df['min'].min():.2f}\n"
                f"Max: {stats_df['max'].max():.2f}\n"
                f"Total Samples: {stats_df['count'].sum()}"
            )

            # Position text box in upper right
            ax.text(
                0.95,
                0.95,
                summary_text,
                transform=ax.transAxes,
                verticalalignment="top",
                horizontalalignment="right",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
            )

            ax.set_title(f"Comparison of {title} Across Plates")
            ax.set_ylabel(signal_column)
            ax.set_xlabel("Plate ID")

            # Rotate x-tick labels if needed
            if len(unique_plates) > 5 or any(len(str(p)) > 10 for p in unique_plates):
                plt.xticks(rotation=45, ha="right")

            # Create a sanitized filename
            safe_signal_column = self.sanitize_filename(signal_column)

            # Save figure
            output_file = os.path.join(
                output_dir, f"plate_comparison_{safe_signal_column}.{save_format}"
            )
            fig.savefig(output_file, dpi=dpi, bbox_inches="tight")
            print(f"Generated plate comparison chart: {output_file}")

            if show_plot:
                plt.show()
            else:
                plt.close(fig)
        except Exception as e:
            error_msg = f"Error generating plate comparison: {str(e)}"
            self.errors.append(error_msg)
            print(error_msg)
            traceback.print_exc()

    def calculate_conversion_ratio(
        self,
        product_column: str = "product-1-EIC-area-Positive",
        reactant_column: str = "reactant-1-EIC-area-Positive",
        ratio_column_name: str = "conversion-ratio",
    ) -> bool:
        """
        Calculate conversion ratio between product and reactant

        Args:
            product_column: Column containing product area values
            reactant_column: Column containing reactant area values
            ratio_column_name: Name for the new ratio column to create

        Returns:
            True if calculation was successful, False otherwise
        """
        if product_column in self.data.columns and reactant_column in self.data.columns:
            try:
                self.data[ratio_column_name] = pd.to_numeric(
                    self.data[product_column], errors="coerce"
                ) / pd.to_numeric(self.data[reactant_column], errors="coerce")
                print(
                    f"Calculated conversion ratio as {product_column} / {reactant_column}"
                )
                return True
            except Exception as e:
                print(f"Error calculating conversion ratio: {str(e)}")
                return False
        else:
            missing = []
            if product_column not in self.data.columns:
                missing.append(product_column)
            if reactant_column not in self.data.columns:
                missing.append(reactant_column)
            print(
                f"Cannot calculate conversion ratio - missing columns: {', '.join(missing)}"
            )
            return False

    def generate_all_visualizations(
        self, metrics: list = None, show_plots: bool = False
    ) -> None:
        """
        Generate all standard visualizations for the given metrics

        Args:
            metrics: List of (column, title, colormap) tuples
            show_plots: Whether to display plots interactively
        """
        if metrics is None:
            metrics = [
                ("product-1-max-EIC-signal-Positive", "Product 1 Signal", "rocket"),
                ("product-1-EIC-area-Positive", "Product 1 Area", "plasma"),
                ("reactant-1-max-EIC-signal-Positive", "Reactant 1 Signal", "inferno"),
                ("reactant-1-EIC-area-Positive", "Reactant 1 Area", "magma"),
            ]

        # Calculate conversion ratio
        if self.calculate_conversion_ratio():
            metrics.append(("conversion-ratio", "Conversion Ratio", "coolwarm"))

        # Generate heatmaps and comparisons
        for column, title, cmap in metrics:
            if column in self.data.columns:
                try:
                    print(f"Generating heatmap for {title}...")
                    self.generate_heatmap(
                        signal_column=column,
                        title=title,
                        colormap=cmap,
                        show_plot=show_plots,
                    )

                    if "plate-ID" in self.data.columns:
                        print(f"Generating plate comparison for {title}...")
                        self.generate_plate_comparison(
                            signal_column=column,
                            title=title,
                            colormap=cmap,
                            show_plot=show_plots,
                        )
                except Exception as e:
                    error_msg = (
                        f"Error generating visualizations for {column}: {str(e)}"
                    )
                    self.errors.append(error_msg)
                    print(error_msg)
                    traceback.print_exc()
            else:
                print(f"Column '{column}' not found in data, skipping visualization")

    @staticmethod
    def sanitize_filename(name: str) -> str:
        """Convert a string to a safe filename"""
        # Replace characters that might cause issues
        replacements = {
            "-": "_",
            " ": "_",
            "+": "plus",
            "/": "_",
            "\\": "_",
            ":": "_",
            "*": "_",
            "?": "_",
            '"': "",
            "<": "_",
            ">": "_",
            "|": "_",
        }

        if not isinstance(name, str):
            name = str(name)

        for char, replacement in replacements.items():
            name = name.replace(char, replacement)

        return name
