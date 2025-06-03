from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import ntpath
import os
import psutil
from logging_config import get_logger

# Get logger for this module
logger = get_logger(__name__)


def get_mol(smiles: str) -> None:
    """
    Creates mol object from target compound smiles
    Args:
        smiles: compound smiles
    """
    return Chem.MolFromSmiles(smiles)


def get_smiles(mol) -> str:
    """
    Creates SMILES string from target compound molß
    Args:
        mol: compound mol
    """
    return Chem.MolToSmiles(mol)


def get_MW(mol) -> float:
    """
    Calculates molecular weight of target compound
    Args:
        mol: mol object target compound
    """
    return round(Descriptors.MolWt(mol))


def get_path_leaf(path):
    """
    Linux and Windows compatible path splitter. Returns
    the final bit at end of path
        Args: path to split
    """
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def create_molecule_svg(mol, filepath=None):
    """
    Creates svg image of rdkit mol and saves to file

    Args:
        mol: RDKit molecule object
        filepath: Path where the SVG should be saved (optional)

    Returns:
        The SVG content as a string
    """
    if mol is None:
        error_msg = "No molecule provided to create_molecule_svg"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Create the SVG
    logger.debug("Generating molecule SVG")
    compound_image = rdMolDraw2D.MolDraw2DSVG(824, 556)
    compound_image.drawOptions().padding = 0
    compound_image.DrawMolecule(mol)
    compound_image.FinishDrawing()
    svg_content = compound_image.GetDrawingText()

    # Save to file if filepath is provided
    if filepath:
        # Make sure the directory exists
        os.makedirs(os.path.dirname(filepath), exist_ok=True)

        with open(filepath, "w") as f:
            f.write(svg_content)
        logger.debug(f"Saved molecule SVG to {filepath}")
    else:
        # Use the old default path for backward compatibility
        default_path = "../tmpimages/molecule.svg"
        os.makedirs(os.path.dirname(default_path), exist_ok=True)

        with open(default_path, "w") as f:
            f.write(svg_content)
        logger.debug(f"Saved molecule SVG to default path {default_path}")

    # Return the SVG content in case it's needed
    return svg_content


def create_report_directory_structure(base_dir: str) -> dict:
    """
    Create a structured directory hierarchy for MS analysis reports

    Args:
        base_dir: Base directory path for reports

    Returns:
        Dictionary with paths to all report subdirectories
    """
    import os

    # Create the base directory if it doesn't exist
    os.makedirs(base_dir, exist_ok=True)

    # Define the directory structure
    dirs = {
        "plate_comparisons": os.path.join(base_dir, "plate-comparisons"),
        "positive": {
            "root": os.path.join(base_dir, "positive-mode"),
            "interactive": os.path.join(base_dir, "positive-mode", "interactive-plots"),
            "static": os.path.join(base_dir, "positive-mode", "static-reports"),
            "heatmaps": os.path.join(base_dir, "positive-mode", "heatmaps"),
        },
        "negative": {
            "root": os.path.join(base_dir, "negative-mode"),
            "interactive": os.path.join(base_dir, "negative-mode", "interactive-plots"),
            "static": os.path.join(base_dir, "negative-mode", "static-reports"),
            "heatmaps": os.path.join(base_dir, "negative-mode", "heatmaps"),
        },
    }

    # Create all directories
    for mode in ["positive", "negative"]:
        for dir_path in dirs[mode].values():
            os.makedirs(dir_path, exist_ok=True)

    os.makedirs(dirs["plate_comparisons"], exist_ok=True)

    return dirs


def monitor_memory(label="Current", logger=None):
    """
    Log current memory usage
    
    Args:
        label: Description of the current monitoring point
        logger: Logger to use (defaults to utils logger if None)
        
    Returns:
        Current memory usage in MB
    """
    if logger is None:
        logger = get_logger(__name__)
        
    try:
        process = psutil.Process()
        mem_info = process.memory_info()
        memory_mb = mem_info.rss / (1024 * 1024)
        logger.info(f"Memory usage ({label}): {memory_mb:.1f} MB")
        
        # Also log system memory if available
        system_memory = psutil.virtual_memory()
        logger.info(f"System memory: {system_memory.percent}% used, {system_memory.available / (1024 * 1024 * 1024):.1f} GB available")
        
        return memory_mb
    except ImportError:
        logger.warning("Memory monitoring requires psutil. Install with 'pip install psutil'")
        return 0
    except Exception as e:
        logger.warning(f"Error monitoring memory: {str(e)}")
        return 0
