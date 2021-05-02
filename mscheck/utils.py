from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from PIL import Image
from io import BytesIO
import matplotlib.pyplot as plt
import ntpath


def get_mol(smiles: str) -> None:
    """
    Creates mol object from target compound smiles
    Args:
        smiles: compound smiles
    """
    return Chem.MolFromSmiles(smiles)


def get_smiles(mol) -> str:
    """
    Creates SMILES string from target compound mol
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
    return round(Chem.Descriptors.MolWt(mol))


def path_leaf(path):
    """
    Linux and Windows compatible path splitter. Returns
    the final bit at end of path
        Args: path to split
    """
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
