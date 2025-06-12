"""UV chromatogram from .mzML files"""

from pyopenms import *
import numpy as np
import re
import logging

logger = logging.getLogger(__name__)


class UVSpectrum(object):
    """
    Creates UV chromatogram object using the pyopenms package
    """

    def __init__(self, mzMLfilepath: str):
        """
        UVSpectrum constructor
        Args:
            mzMLfilepath (str): path to .mzML file
        """
        self._filepath = mzMLfilepath
        self._exp = MSExperiment()
        MzMLFile().load(self._filepath, self._exp)
        
        # Extract UV data
        self.has_uv = False
        self.chromatograms = []
        self.wavelengths = []
        
        # Try to extract UV chromatograms
        self._extract_uv_chromatograms()
        
    def _extract_uv_chromatograms(self):
        """
        Extract UV chromatograms from the mzML file
        """
        try:
            # Track whether we found UV chromatograms
            found_uv = False
            
            # Get all chromatograms
            for i in range(self._exp.getNrChromatograms()):
                chrom = self._exp.getChromatogram(i)
                
                # Try to determine if this is a UV chromatogram
                chrom_id = chrom.getNativeID()
                logger.debug(f"Processing chromatogram ID: {chrom_id}")
                
                # Check if this looks like a UV chromatogram by ID
                is_uv = False
                wavelength = None
                
                if any(term in chrom_id.upper() for term in ["UV", "DAD", "PDA", "WAVELENGTH"]):
                    is_uv = True
                    # Use the helper function to extract wavelength
                    wavelength = self._extract_wavelength_from_id(chrom_id)
                
                # If identified as UV, add to our collection
                if is_uv:
                    found_uv = True
                    self.chromatograms.append(chrom)
                    self.wavelengths.append(wavelength if wavelength else 0)
                    logger.debug(f"Added UV chromatogram with wavelength: {wavelength}")
        
            # If we found UV chromatograms, set flag
            if found_uv:
                self.has_uv = True
                logger.info(f"Found {len(self.chromatograms)} UV chromatograms")
                
                # Log extracted wavelengths if available
                if any(w > 0 for w in self.wavelengths):
                    valid_wavelengths = [w for w in self.wavelengths if w > 0]
                    logger.info(f"Extracted wavelengths: {valid_wavelengths}")
            else:
                logger.info("No UV chromatograms found in the file")
                
        except Exception as e:
            logger.error(f"Error extracting UV chromatograms: {str(e)}")
            
    def get_chromatogram(self, target_wavelength=None):
        """
        Get UV chromatogram at or closest to the specified wavelength
        
        Args:
            target_wavelength: Target wavelength in nm, or None for the first chromatogram
            
        Returns:
            tuple: (retention_times, intensities) or (None, None) if no chromatogram found
        """
        if not self.has_uv or not self.chromatograms:
            return None, None
            
        # If no wavelength specified, return the first chromatogram
        if target_wavelength is None:
            chrom = self.chromatograms[0]
            rt, intensity = chrom.get_peaks()
            return rt / 60.0, intensity  # Convert to minutes
        
        # Find chromatogram closest to the target wavelength
        best_idx = 0
        best_diff = float('inf')
        
        for i, wavelength in enumerate(self.wavelengths):
            if wavelength > 0:  # Skip unknown wavelengths
                diff = abs(wavelength - target_wavelength)
                if diff < best_diff:
                    best_diff = diff
                    best_idx = i
        
        # Get the selected chromatogram
        chrom = self.chromatograms[best_idx]
        rt, intensity = chrom.get_peaks()
        
        return rt / 60.0, intensity  # Convert to minutes
        
    def get_wavelengths(self):
        """
        Get the available wavelengths for UV chromatograms
        
        Returns:
            list: List of available wavelengths (0 means wavelength not known)
        """
        return self.wavelengths
        
    def get_multi_wavelength_data(self):
        """
        Get data for all available wavelengths
        
        Returns:
            dict: Dictionary with wavelengths, retention times and intensities
        """
        result = {
            'wavelengths': [],
            'retention_times': [],
            'intensities': []
        }
        
        for i, chrom in enumerate(self.chromatograms):
            rt, intensity = chrom.get_peaks()
            result['wavelengths'].append(self.wavelengths[i])
            result['retention_times'].append(rt / 60.0)  # Convert to minutes
            result['intensities'].append(intensity)
            
        return result
        
    def get_max_absorption_chromatogram(self):
        """
        Generate a chromatogram of maximum absorption across all wavelengths
        
        Returns:
            tuple: (retention_times, max_intensities) or (None, None) if no data
        """
        if not self.has_uv or not self.chromatograms:
            return None, None
            
        # Get the first chromatogram to establish retention time points
        rt_base, _ = self.chromatograms[0].get_peaks()
        rt = rt_base / 60.0  # Convert to minutes
        
        # Check if all chromatograms have the same retention time points
        same_rt = True
        for chrom in self.chromatograms[1:]:
            rt_points, _ = chrom.get_peaks()
            if len(rt_points) != len(rt_base) or not np.allclose(rt_points, rt_base):
                same_rt = False
                break
        
        if same_rt and len(self.chromatograms) > 1:
            # Create matrix of all intensities
            intensity_matrix = np.vstack([chrom.get_peaks()[1] for chrom in self.chromatograms])
            
            # Get maximum intensity at each retention time point
            max_intensity = np.max(intensity_matrix, axis=0)
            return rt, max_intensity
        elif len(self.chromatograms) > 0:
            # If retention times differ, just return the first chromatogram
            _, intensity = self.chromatograms[0].get_peaks()
            return rt, intensity
        else:
            return None, None
        

    def _extract_wavelength_from_id(self, chrom_id):
        """
        Extract wavelength from chromatogram ID using multiple approaches
        
        Args:
            chrom_id (str): Chromatogram ID string
            
        Returns:
            float or None: Extracted wavelength if found, None otherwise
        """
        wavelength = None
        
        # Method 1: Look for common patterns with regex
        patterns = [
            r'Sig=(\d+)',           # Sig=220 format
            r'WAVELENGTH[=:](\d+)', # WAVELENGTH=220 or WAVELENGTH:220
            r'WVL[=:](\d+)',        # WVL=220
            r'(\d+)nm',             # 220nm format
            r'[\s_](\d{3})[\s_]'    # Space/underscore + 3 digits + space/underscore
        ]
        
        for pattern in patterns:
            match = re.search(pattern, chrom_id, re.IGNORECASE)
            if match:
                try:
                    candidate = float(match.group(1))
                    if 190 <= candidate <= 800:  # Typical UV range
                        logger.debug(f"Found wavelength {candidate} from pattern {pattern} in {chrom_id}")
                        return candidate
                except ValueError:
                    pass
                    
        # Method 2: Original approach - look for standalone numeric parts
        for part in chrom_id.split():
            if part.replace('.', '', 1).isdigit():
                try:
                    candidate = float(part)
                    if 190 <= candidate <= 800:  # Typical UV range
                        logger.debug(f"Found wavelength {candidate} from split part in {chrom_id}")
                        return candidate
                except ValueError:
                    pass
        
        return wavelength
        

# import matplotlib.pyplot as plt

# # Create UV spectrum object
# uv = UVSpectrum(mzMLfilepath="/Users/bvh64415/Library/CloudStorage/OneDrive-DiamondLightSourceLtd/MSCheck-dev/Human conversion estimates/products MS data/open_source_with_uv/PRODUCT-2.mzML")

# # Check if UV data was found
# if uv.has_uv:
#     # Available wavelengths
#     print(f"Available wavelengths: {uv.get_wavelengths()}")
    
#     # Get chromatogram at 254nm (or closest available)
#     wavelength = 280
#     rt, intensity = uv.get_chromatogram(target_wavelength=wavelength)
    
#     # Plot the chromatogram
#     plt.figure(figsize=(10, 6))
#     plt.plot(rt, intensity)
#     plt.xlabel("Retention Time (min)")
#     plt.ylabel("Absorbance")
#     plt.title(f"UV Chromatogram ({wavelength} nm)")
#     plt.grid(True)
#     plt.show()
    
#     # Get maximum absorption chromatogram across all wavelengths
#     rt_max, intensity_max = uv.get_max_absorption_chromatogram()
    
#     # Plot maximum absorption
#     plt.figure(figsize=(10, 6))
#     plt.plot(rt_max, intensity_max)
#     plt.xlabel("Retention Time (min)")
#     plt.ylabel("Maximum Absorbance")
#     plt.title("Maximum UV Absorption Across All Wavelengths")
#     plt.grid(True)
#     plt.show()
# else:
#     print("No UV data found in the mzML file")