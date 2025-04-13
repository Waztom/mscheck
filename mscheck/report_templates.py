"""HTML and JavaScript templates for MS report generation"""

import os
import json

def create_interactive_report_html(report_title, plot_figure, mass_spectra, 
                                  json_encoder=None, available_reports=None, 
                                  current_report_index=0):
    """
    Create a complete interactive HTML report using external template files
    
    Args:
        report_title: Title for the report
        plot_figure: Plotly figure object to display
        mass_spectra: List of mass spectra data
        json_encoder: Optional custom JSON encoder for NumPy arrays
        available_reports: List of all available reports for navigation
        current_report_index: Index of current report in the available_reports list
        
    Returns:
        Complete HTML content as string
    """
    # Convert mass spectra to JSON
    mass_spectra_json = json.dumps(mass_spectra, cls=json_encoder)
    
    # Extract figure data properly
    data = []
    for trace in plot_figure.data:
        data.append(trace.to_plotly_json())
        
    layout = plot_figure.layout.to_plotly_json()
    
    # Create a clean JSON string with the extracted data using NumpyEncoder
    plot_json = json.dumps({"data": data, "layout": layout}, cls=json_encoder)
    
    # Handle reports navigation data
    if available_reports is None:
        # Default to just the current report if no list provided
        available_reports = [{
            "title": report_title,
            "path": "#"
        }]
    
    # Convert navigation data to JSON
    reports_json = json.dumps(available_reports, cls=json_encoder)
    
    # Read the HTML template
    template_path = os.path.join(os.path.dirname(__file__), 'templates', 'report_template.html')
    with open(template_path, 'r') as f:
        template_content = f.read()
    
    # Read the spectrum lookup script
    script_path = os.path.join(os.path.dirname(__file__), 'templates', 'spectrum_lookup.js')
    with open(script_path, 'r') as f:
        spectrum_script = f.read()
    
    # Replace placeholders
    html_content = template_content.replace('{{REPORT_TITLE}}', report_title)
    html_content = html_content.replace('{{PLOT_JSON}}', plot_json)
    html_content = html_content.replace('{{MASS_SPECTRA_JSON}}', mass_spectra_json)
    html_content = html_content.replace('{{AVAILABLE_REPORTS}}', reports_json)
    html_content = html_content.replace('{{CURRENT_REPORT_INDEX}}', str(current_report_index))
    html_content = html_content.replace('{{SPECTRUM_SCRIPT}}', spectrum_script)
    
    return html_content