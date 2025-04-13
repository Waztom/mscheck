"""Generate report function"""

# Future imports
from __future__ import annotations

# Standard library imports
import os
import shutil
from pathlib import Path
from typing import Dict, List

# Third-party imports
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from svgutils.compose import Figure, SVG

# Local imports
from logging_config import get_logger
from utils import create_molecule_svg

# Get logger for this module
logger = get_logger(__name__)


class MSReport:
    """Generates detailed reports for MS analysis results"""

    # Define default plot parameters
    DEFAULT_PARAMS = {
        "font.weight": "bold",
        "legend.fontsize": "small",
        "figure.figsize": (12, 12),
        "axes.labelsize": "large",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "axes.titlesize": "x-large",
        "xtick.labelsize": "large",
        "ytick.labelsize": "large",
        "figure.titleweight": "bold",
        "figure.subplot.bottom": 0.11,
        "figure.subplot.hspace": 0.2,
        "figure.subplot.left": 0.125,
        "figure.subplot.right": 0.9,
        "figure.subplot.top": 0.88,
        "figure.subplot.wspace": 0.2,
        "figure.constrained_layout.use": True,
    }

    def __init__(self, output_dir="reports", temp_dir=None):
        """
        Initialize report generator
        """
        from logging_config import get_logger
        
        # Use the main mscheck logger, not a class-specific one
        self.logger = get_logger("mscheck")  # This will use the same logger as BulkAnalyser
        
        self.output_dir = output_dir
        self.temp_dir = temp_dir or os.path.join(output_dir, "temp")
        self.errors = []
        
        # Create directories
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.temp_dir, exist_ok=True)
        
        self.logger.info(f"MSReport initialized with output directory: {self.output_dir}")

    def __del__(self):
        """Clean up temporary files when instance is deleted"""
        try:
            if Path(self.temp_dir).exists():
                shutil.rmtree(self.temp_dir)
                self.logger.debug(f"Cleaned up temporary directory: {self.temp_dir}")
        except Exception as e:
            self.logger.warning(f"Failed to clean up temporary directory: {str(e)}")
            pass

    def create_compound_report(
        self,
        msmode: str,
        RT_values: list,
        TIC_values: list,
        compound_name: str,
        mol,  # RDKit molecule
        analysedata: dict = None,
    ) -> str:
        """
        Create a detailed report for a single compound

        Args:
            msmode: "Positive" or "Negative"
            RT_values: List of retention time values
            TIC_values: List of total ion count values
            compound_name: Name of the compound
            mol: RDKit molecule object
            analysedata: Dictionary containing analysis results

        Returns:
            Path to the generated report file
        """
        try:
            self.logger.info(f"Creating compound report for: {compound_name}")
            
            # Determine ion mode
            ion_mode = "+" if msmode == "Positive" else "-"

            # Determine if we have matching ions
            no_plots = 0 if analysedata is None else len(analysedata.get("ions", []))
            self.logger.debug(f"Number of ion plots to create: {no_plots}")

            if no_plots == 0:
                # Create simple report if no ions found
                self.logger.info(f"No ions found for {compound_name}, creating simple report")
                output_path = self._create_simple_report(
                    RT_values, TIC_values, compound_name, ion_mode, mol
                )
            else:
                # Create detailed report with ion matches
                self.logger.info(f"Found {no_plots} ions for {compound_name}, creating detailed report")
                output_path = self._create_detailed_report(
                    RT_values,
                    TIC_values,
                    compound_name,
                    ion_mode,
                    mol,
                    analysedata,
                    no_plots,
                )

            self.logger.info(f"Successfully created report at: {output_path}")
            return output_path

        except Exception as e:
            error_msg = f"Error creating report for {compound_name}: {str(e)}"
            self.errors.append(error_msg)
            self.logger.error(error_msg)
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _create_simple_report(self, RT_values, TIC_values, compound_name, ion_mode, mol):
        """Create a simple report when no masses are found"""
        self.logger.debug(f"Creating simple report for {compound_name}")
        
        fig, ax = plt.subplots(figsize=(12, 4))
        fig.suptitle(f"MSCheck report: Mass not found for {compound_name}", size=16)

        ax.plot(RT_values, TIC_values, color="k", zorder=-1)
        ax.set_xlabel("Retention time (min)")
        ax.set_ylabel("Total ion count (TIC)")
        ax.set_xlim([-0.9, ax.get_xlim()[1]])
        ax.set_ylim([0, ax.get_ylim()[1]])

        # Save temporary files
        plot_path = f"{self.temp_dir}/plot.svg"
        fig.savefig(plot_path, transparent=True)
        plt.close(fig)
        self.logger.debug(f"Saved plot to {plot_path}")

        # Create molecule SVG
        mol_path = f"{self.temp_dir}/molecule.svg"
        try:
            create_molecule_svg(mol, filepath=mol_path)
            self.logger.debug(f"Created molecule SVG at {mol_path}")
        except TypeError:
            # If filepath parameter doesn't work, try output_path or other parameter name
            self.logger.warning("Parameter name issue with create_molecule_svg, trying alternate parameter name")
            try:
                create_molecule_svg(mol, output_path=mol_path)
                self.logger.debug(f"Created molecule SVG at {mol_path} with alternate parameter")
            except Exception as e:
                self.logger.error(f"Failed to create molecule SVG: {str(e)}")
                raise

        # Combine into final report
        output_path = f"{self.output_dir}/{compound_name}-report.svg"
        Figure(
            "29cm",
            "40cm",
            SVG(mol_path).scale(0.004).move(3, 3),
            SVG(plot_path).scale(0.03),
        ).save(output_path)
        self.logger.debug(f"Saved final report to {output_path}")

        return output_path

    def _create_detailed_report(
        self,
        RT_values,
        TIC_values,
        compound_name,
        ion_mode,
        mol,
        analysedata,
        no_plots,
    ):
        """Create a detailed report with ion matches"""
        self.logger.debug(f"Creating detailed report for {compound_name}")
        
        # Create figure with subplots
        fig, ax = plt.subplots((no_plots * 2) + 1, figsize=(12, no_plots * 4 + 4))
        fig.tight_layout(pad=4.2)
        fig.suptitle(f"MSCheck report: Mass found for {compound_name}", size=16)
        fig.align_ylabels()

        # Plot TIC
        ax[0].plot(RT_values, TIC_values, color="k", zorder=-1)
        ax[0].set_xlabel("Retention time (min)")
        ax[0].set_ylabel("Total ion count (TIC)")
        ax[0].set_xlim([-0.9, ax[0].get_xlim()[1]])
        ax[0].set_ylim([0, ax[0].get_ylim()[1]])

        subplot = 1
        colors = iter(cm.rainbow(np.linspace(0, 1, no_plots * 2)))

        # Add ion data plots
        for (
            EIC_data,
            max_mz_match,
            ion_found,
            RT_match_values,
            TIC_match_values,
            mz_strongest,
        ) in zip(
            analysedata["EIC_data"],
            analysedata["max_mz_match"],
            analysedata["ions"],
            analysedata["RT"],
            analysedata["TIC"],
            analysedata["mz_strongest"],
        ):
            max_match_label = "Dominant" if max_mz_match else "Not dominant"
            ion_name = ion_found[0].strip("[]")
            mz_masses_max, mz_intensities_max, max_index = mz_strongest

            color_matches = next(colors)

            # Plot non-max points
            ax[0].scatter(
                [RT for i, RT in enumerate(RT_match_values) if i != max_index],
                [TIC for i, TIC in enumerate(TIC_match_values) if i != max_index],
                color=color_matches,
                s=45.0,
                linewidth=3,
                marker="x",
                zorder=1,
                label=f"{max_match_label} mz match for M{ion_mode}{ion_name}",
            )
            ax[0].legend(loc="upper right")

            # Plot max point
            RT_max = RT_match_values[max_index]
            TIC_max = TIC_match_values[max_index]
            color_max = next(colors)
            ax[0].scatter(
                RT_max,
                TIC_max,
                color=color_max,
                s=45.0,
                linewidth=3,
                marker="o",
                zorder=1,
                label=f"Strongest mz match pattern for M{ion_mode}{ion_name}",
            )
            ax[0].legend(loc="upper right")

            # Create mass spectrum plot
            markerline, stemline, baseline = ax[subplot].stem(
                mz_masses_max, mz_intensities_max
            )
            plt.setp(markerline, "markerfacecolor", color_max, "markersize", 10)
            plt.setp(stemline, "color", "k")
            plt.setp(baseline, "color", "k")

            ax[subplot].set_title(
                f"Stongest mz pattern matching M{ion_mode}{ion_name} "
                f"({ion_found[1]}) at RT: {np.round(RT_max, 1)} (min)"
            )
            ax[subplot].set_xlabel("m/z (Da)")
            ax[subplot].set_ylabel("Relative intensity")
            ax[subplot].set_xlim(
                [ax[subplot].get_xlim()[0], ax[subplot].get_xlim()[1] + 10]
            )

            # Add m/z annotations
            for i, j in zip(mz_masses_max, mz_intensities_max):
                ax[subplot].annotate(
                    str(i),
                    xy=(i, j),
                    textcoords="offset points",
                    xytext=(5.5, -3.5),
                    ha="left",
                )

            subplot += 1

            # Create EIC plot
            ax[subplot].set_title(f"EIC for M{ion_mode}{ion_name}")
            ax[subplot].set_xlabel("Retention time (min)")
            ax[subplot].set_ylabel("Extracted ion intensity")
            ax[subplot].set_xlim([-0.9, ax[0].get_xlim()[1]])

            ax[subplot].scatter(
                [data[0] for data in EIC_data],
                [data[1] for data in EIC_data],
                color=color_max,
            )

            subplot += 1

        # Save temporary files
        plot_path = f"{self.temp_dir}/plot.svg"
        fig.savefig(plot_path, transparent=True)
        plt.close(fig)
        self.logger.debug(f"Saved plot to {plot_path}")

        # Create molecule SVG
        mol_path = f"{self.temp_dir}/molecule.svg"
        try:
            create_molecule_svg(mol, filepath=mol_path)
            self.logger.debug(f"Created molecule SVG at {mol_path}")
        except TypeError:
            self.logger.warning("Parameter name issue with create_molecule_svg, trying alternate parameter name")
            try:
                create_molecule_svg(mol, output_path=mol_path)
                self.logger.debug(f"Created molecule SVG at {mol_path} with alternate parameter")
            except Exception as e:
                self.logger.error(f"Failed to create molecule SVG: {str(e)}")
                raise

        # Combine into final report
        output_path = f"{self.output_dir}/{compound_name}-report.svg"
        Figure(
            "29cm",
            "40cm",
            SVG(mol_path).scale(0.0035).move(3.4, 2.0),
            SVG(plot_path).scale(0.03),
        ).save(output_path)
        self.logger.debug(f"Saved final report to {output_path}")

        return output_path

    def create_annotated_tic_report(
        self,
        RT_values: list,
        TIC_values: list,
        compounds: List[Dict],
        report_title: str = "Annotated TIC Report",
        html_output: bool = True,
        svg_layout: str = "standard",
        available_reports: List[Dict] = None,
        current_report_index: int = 0
    ) -> str:
        """
        Create an annotated TIC plot with molecule structures and mass spectra
        
        Args:
            RT_values: Retention time values for the master TIC
            TIC_values: TIC intensity values
            compounds: List of dictionaries with compound data
            report_title: Title for the report
            html_output: Whether to generate HTML (True) or SVG (False)
            svg_layout: Layout for SVG output - "standard" or "triple"
            available_reports: List of all available reports for navigation
            current_report_index: Index of current report in available_reports
            
        Returns:
            Path to the generated report file
        """
        try:
            # For HTML output
            if html_output:
                self.logger.info("Creating interactive HTML report")
                return self._create_interactive_report(
                    RT_values=RT_values,
                    TIC_values=TIC_values,
                    compounds=compounds,
                    report_title=report_title,
                    available_reports=available_reports,
                    current_report_index=current_report_index
                )
            
            # For SVG output
            elif svg_layout == "triple":
                return self._create_triple_section_svg(RT_values, TIC_values, compounds, report_title)
            
            # For standard SVG output
            else:
                return self._create_standard_svg(RT_values, TIC_values, compounds, report_title)
                
        except Exception as e:
            error_msg = f"Error creating annotated TIC report: {str(e)}"
            self.errors.append(error_msg)
            self.logger.error(error_msg)
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _create_interactive_report(
        self,
        RT_values: list,
        TIC_values: list,
        compounds: List[Dict],
        report_title: str,
        available_reports: List[Dict] = None,
        current_report_index: int = 0
    ) -> str:
        """Create interactive HTML report with complete MZ data for all retention times"""
        import plotly.graph_objects as go
        import numpy as np
        from report_templates import create_interactive_report_html
        from numpyencoder import NumpyEncoder
        
        # Log the original data size
        self.logger.info(f"Processing {len(RT_values)} data points")
        
        # Convert to numpy arrays and ensure proper types
        RT_values = np.array(RT_values, dtype=float)
        TIC_values = np.array(TIC_values, dtype=float)
        

        # Create figure
        fig = go.Figure()
        
        # 1. Add visible TIC line trace first (underneath everything)
        fig.add_trace(go.Scatter(
            x=RT_values,
            y=TIC_values,
            mode='lines',
            name='TIC',
            line=dict(color='black', width=1.5),
            hoverinfo='skip'
        ))
        
        # 3. Add interactive points with clear visibility
        fig.add_trace(go.Scatter(
            x=RT_values,
            y=TIC_values,
            mode='markers',
            name='TIC Points',
            marker=dict(
                size=6,
                color='rgba(65, 105, 225, 0.7)',  # Royal blue with transparency
                line=dict(width=1, color='rgba(25, 25, 112, 0.5)')  # Subtle border
            ),
            hovertemplate="RT: %{x:.2f} min<br>TIC: %{y:,.0f}<extra></extra>",
            customdata=[[rt, -1] for rt in RT_values],  # -1 indicates not a compound
            showlegend=False
        ))
        
        # 5. Prepare mass spectra data - COMBINING ALL MZ DATA SOURCES
        mass_spectra = []
        
        # First look for all_mz_data in the compounds
        all_mz_data_found = False
        for compound in compounds:
            if 'all_mz_data' in compound and compound['all_mz_data'] is not None:
                all_mz_data_found = True
                all_mz_data = compound['all_mz_data']
                self.logger.info(f"Using all_mz_data from compound {compound.get('compound_name')}")
                break
        
        # If found, add MZ data for all RT points
        if all_mz_data_found and all_mz_data:
            self.logger.info(f"Processing {len(all_mz_data)} MZ spectra for interactive plot")
            
            for i, (rt, mz_spectrum) in enumerate(zip(RT_values, all_mz_data)):
                if isinstance(mz_spectrum, tuple) and len(mz_spectrum) >= 2:
                    mz_values, intensity_values = mz_spectrum
                    
                    if (isinstance(mz_values, (list, np.ndarray)) and 
                        isinstance(intensity_values, (list, np.ndarray)) and
                        len(mz_values) > 0 and len(intensity_values) > 0):
                        
                        # Add MZ spectrum for this RT point
                        mass_spectra.append({
                            "index": -i-1,  # Use negative indices for scan points
                            "name": f"Scan at RT={rt:.2f}",
                            "compound_type": "scan",
                            "rt": rt,
                            "mz": mz_values.tolist() if hasattr(mz_values, 'tolist') else list(mz_values),
                            "intensity": intensity_values.tolist() if hasattr(intensity_values, 'tolist') else list(intensity_values),
                            "color": "gray"
                        })
        
        # Then add strongest matches from compound data (with colors)
        for idx, compound in enumerate(compounds):
            name = compound.get("compound_name", "Unknown")
            compound_type = compound.get("compound_type", "unknown")
            rt_max = compound.get("rt_max")
            mz_data = compound.get("mz_strongest")
            
            if isinstance(mz_data, tuple) and len(mz_data) >= 2:
                mz_values, intensity_values = mz_data
                
                if (isinstance(mz_values, (list, np.ndarray)) and 
                    isinstance(intensity_values, (list, np.ndarray)) and
                    len(mz_values) > 0 and len(intensity_values) > 0):
                    
                    # Define color based on compound type
                    color = {
                        "reactant": "red",
                        "product": "green", 
                        "internal-std": "blue",
                        "intermediate": "purple"
                    }.get(compound_type, "orange")
                    
                    # Add to mass spectra data
                    mass_spectra.append({
                        "index": idx,
                        "name": name,
                        "compound_type": compound_type,
                        "rt": rt_max,
                        "mz": mz_values.tolist() if hasattr(mz_values, 'tolist') else list(mz_values),
                        "intensity": intensity_values.tolist() if hasattr(intensity_values, 'tolist') else list(intensity_values),
                        "color": color
                    })
                    
                    # Find closest index in TIC data
                    closest_idx = min(range(len(RT_values)), key=lambda i: abs(RT_values[i] - rt_max))
                    
                    # Add compound marker
                    fig.add_trace(go.Scatter(
                        x=[rt_max],
                        y=[TIC_values[closest_idx]],
                        mode='markers',
                        marker=dict(
                            color=color,
                            size=12,
                            symbol='diamond',
                            line=dict(color='black', width=1.5)
                        ),
                        name=name,
                        hovertext=f"<b>{name}</b><br>RT: {rt_max:.2f}<br>Type: {compound_type}",
                        hoverinfo="text",
                        customdata=[[rt_max, idx]],
                        showlegend=True
                    ))
        
        # 6. Configure layout with clearer appearance
        fig.update_layout(
            title={
                'text': report_title,
                'y': 0.95,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            xaxis_title="Retention Time (min)",
            yaxis_title="Total Ion Count (TIC)",
            height=600,
            hovermode="closest",
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.15,
                xanchor="center",
                x=0.5
            ),
            plot_bgcolor='rgba(240,240,240,0.5)'  # Light gray background for better contrast
        )
        
        # 7. Create HTML report
        html_path = f"{self.output_dir}/{self.sanitize_filename(report_title)}.html"
        
        if mass_spectra:
            # Create complete HTML with navigation
            html_content = create_interactive_report_html(
                report_title=report_title,
                plot_figure=fig,
                mass_spectra=mass_spectra,
                json_encoder=NumpyEncoder,
                available_reports=available_reports,
                current_report_index=current_report_index
            )
            
            # Write HTML file
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_content)
        else:
            # If no mass spectra available, just write the basic figure to HTML
            fig.write_html(html_path)
        
        self.logger.info(f"Saved interactive HTML report to {html_path}")
        return html_path

    def sanitize_filename(self, name: str) -> str:
        """Convert a string to a safe filename"""
        # Replace characters that might cause issues
        replacements = {
            '-': '_',
            ' ': '_',
            '+': 'plus',
            '/': '_',
            '\\': '_',
            ':': '_',
            '*': '_',
            '?': '_',
            '"': '',
            '<': '_',
            '>': '_',
            '|': '_'
        }
        
        if not isinstance(name, str):
            name = str(name)
            
        for char, replacement in replacements.items():
            name = name.replace(char, replacement)
        
        return name


def get_marker_props(signal):
    """Get marker properties based on signal strength"""
    if isinstance(signal, (int, float)) and signal > 50000:
        return {"symbol": "star", "size": 15, "line_width": 2}
    elif isinstance(signal, (int, float)) and signal > 10000:
        return {"symbol": "diamond", "size": 12, "line_width": 1.5}
    else:
        return {"symbol": "circle", "size": 10, "line_width": 1}
