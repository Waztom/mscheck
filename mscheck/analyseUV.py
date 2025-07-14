"""Analyse UV spectrum class"""

from __future__ import annotations
import numpy as np
import logging
from scipy.signal import find_peaks
from uvspectrum import UVSpectrum

class AnalyseUV(UVSpectrum):
    """
    Analyses the UVSpectrum class object to extract peak information,
    calculate areas, and determine optimal wavelengths
    """

    def __init__(self, mzMLfilepath: str):
        """
        Analyse UV spectrum constructor
        
        Args:
            mzMLfilepath (str): path to .mzML file
        """
        super().__init__(mzMLfilepath)
        self.logger = logging.getLogger("AnalyseUV")
        self.logger.info(f"Initialized AnalyseUV for {mzMLfilepath}")
        
    def find_peak_at_rt(self, retention_time: float, tolerance: float = 0.5) -> dict:
        """
        Find UV peaks at the specified retention time within the given tolerance
        
        Args:
            retention_time (float): Target retention time in minutes
            tolerance (float): Time tolerance window in minutes
            
        Returns:
            dict: Peak information including RT, intensity values, and wavelengths
        """
        if not self.has_uv:
            self.logger.warning("No UV chromatograms available for analysis")
            return {'peaks_found': False}
            
        self.logger.info(f"Searching for peak at RT {retention_time} min ± {tolerance} min")
        
        # Define retention time window
        rt_min = retention_time - tolerance
        rt_max = retention_time + tolerance
        
        # Get all multi-wavelength data from UVSpectrum class
        uv_data = self.get_multi_wavelength_data()
        
        # Store data for each wavelength within the RT window
        peak_data = {
            'wavelengths': [],
            'chromatogram_indices': [],
            'rt_indices': [],
            'rt_values': [],
            'intensity_values': [],
            'peaks_found': False
        }
        
        # Process each wavelength
        for i, (wavelength, rt_array, intensity_array) in enumerate(
                zip(uv_data['wavelengths'], uv_data['retention_times'], uv_data['intensities'])):
            
            # Find indices within retention time window
            window_indices = np.where((rt_array >= rt_min) & (rt_array <= rt_max))[0]
            
            if len(window_indices) > 0:
                # Extract RT and intensity values
                rt_values = rt_array[window_indices]
                intensity_values = intensity_array[window_indices]
                
                # Only add if there's actual signal
                if np.max(intensity_values) > 0:
                    peak_data['wavelengths'].append(wavelength)
                    peak_data['chromatogram_indices'].append(i)
                    peak_data['rt_indices'].append(window_indices)
                    peak_data['rt_values'].append(rt_values)
                    peak_data['intensity_values'].append(intensity_values)
                    peak_data['peaks_found'] = True
                    self.logger.debug(f"Found signal at wavelength {wavelength} nm")
                    
        if peak_data['peaks_found']:
            self.logger.info(f"Found peak data at RT {retention_time} min across {len(peak_data['wavelengths'])} wavelengths")
        else:
            self.logger.warning(f"No peak found at RT {retention_time} min within tolerance {tolerance} min")
            
        return peak_data
    
    def calculate_peak_area(self, peak_data: dict, wavelength: float = None) -> dict:
        """
        Calculate the area of a UV peak using the trapezoidal rule
        
        Args:
            peak_data (dict): Peak information from find_peak_at_rt
            wavelength (float, optional): Specific wavelength to calculate area for
                                        If None, calculates area for all wavelengths
        
        Returns:
            dict: Dictionary mapping wavelengths to their respective areas
        """
        if not peak_data.get('peaks_found', False):
            self.logger.warning("No peak data available for area calculation")
            return {}
            
        areas = {}
        
        # Process data for each wavelength
        for i, wl in enumerate(peak_data['wavelengths']):
            # Skip if not the requested wavelength
            if wavelength is not None and wl != wavelength:
                continue
                
            rt_values = peak_data['rt_values'][i]
            intensity_values = peak_data['intensity_values'][i]
            
            # Calculate area using trapezoidal rule
            area = 0.0
            for j in range(1, len(rt_values)):
                # Width of the interval (difference in retention time)
                rt_diff = rt_values[j] - rt_values[j-1]
                
                # Average height (average intensity)
                height_avg = (intensity_values[j] + intensity_values[j-1]) / 2
                
                # Area of the trapezoid
                area += rt_diff * height_avg
                
            areas[wl] = area
            self.logger.debug(f"Area for wavelength {wl} nm: {area:.2f}")
            
        self.logger.info(f"Calculated areas for {len(areas)} wavelengths")
        return areas
    
    def find_optimal_wavelength(self, peak_data: dict) -> tuple:
        """
        Determine which wavelength gives the strongest response (area) for a peak
        
        Args:
            peak_data (dict): Peak information from find_peak_at_rt
            
        Returns:
            tuple: (optimal_wavelength, max_area)
        """
        if not peak_data.get('peaks_found', False):
            self.logger.warning("No peak data available to find optimal wavelength")
            return (None, 0)
            
        # Calculate area for all wavelengths
        areas = self.calculate_peak_area(peak_data)
        
        if not areas:
            return (None, 0)
            
        # Find wavelength with maximum area
        optimal_wavelength = max(areas, key=areas.get)
        max_area = areas[optimal_wavelength]
        
        self.logger.info(f"Optimal wavelength: {optimal_wavelength} nm with area: {max_area:.2f}")
        return (optimal_wavelength, max_area)
    
    def analyse(self, retention_time: float, tolerance: float = 0.5) -> dict:
        """
        Comprehensive analysis of UV data at a specific retention time
        
        Args:
            retention_time (float): Target retention time in minutes
            tolerance (float): Time tolerance window in minutes
            
        Returns:
            dict: Complete analysis results
        """
        # Find peak at specified retention time
        peak_data = self.find_peak_at_rt(retention_time, tolerance)
        
        if not peak_data.get('peaks_found', False):
            self.logger.warning(f"No UV peaks found at RT {retention_time} min ± {tolerance} min")
            return {'found': False}
        
        # Calculate areas for all wavelengths
        areas = self.calculate_peak_area(peak_data)
        
        # Find optimal wavelength
        optimal_wavelength, max_area = self.find_optimal_wavelength(peak_data)
        
        # Find peak information for optimal wavelength
        idx = peak_data['wavelengths'].index(optimal_wavelength)
        rt_values = peak_data['rt_values'][idx]
        intensity_values = peak_data['intensity_values'][idx]
        max_intensity = np.max(intensity_values)
        max_intensity_rt = rt_values[np.argmax(intensity_values)]
        
        analysis_results = {
            'found': True,
            'areas': areas,
            'optimal_wavelength': optimal_wavelength,
            'max_area': max_area,
            'max_intensity': max_intensity,
            'max_intensity_rt': max_intensity_rt,
            'rt_range': (np.min(rt_values), np.max(rt_values)),
            'wavelengths': peak_data['wavelengths'],
        }
        
        self.logger.info(f"Analysis complete. Optimal wavelength: {optimal_wavelength} nm, Max area: {max_area:.2f}")
        return analysis_results
    
# # Initialize analyzer with mzML file
# uv_analyzer = AnalyseUV("/Users/bvh64415/Library/CloudStorage/OneDrive-DiamondLightSourceLtd/FFF-projects/DENV-NS2B3-NS3(MedChemica)/CAR/flavi-t3c-i2a/QC/flavi-t3c-i2a xp00-xp02 with IS/open_source_with_uv/lp02/c0800.mzML")

# # Analyze UV data at retention time 2.5 minutes with 0.2 minute tolerance
# results = uv_analyzer.analyse(retention_time=1.34, tolerance=0.1)

# if results['found']:
#     print(f"Optimal wavelength: {results['optimal_wavelength']} nm")
#     print(f"Maximum area: {results['max_area']}")
#     print(f"Maximum intensity: {results['max_intensity']} at RT {results['max_intensity_rt']} min")
# else:
#     print("No peaks found in the specified retention time window")