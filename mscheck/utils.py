from __future__ import annotations
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
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


def get_path_leaf(path):
    """
    Linux and Windows compatible path splitter. Returns
    the final bit at end of path
        Args: path to split
    """
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def create_molecule_svg(mol: rdkitmol):
    """
    Creates svg image of rdkit mol and saves to file
    """
    compound_image = rdMolDraw2D.MolDraw2DSVG(824, 556)
    compound_image.drawOptions().padding = 0
    compound_image.DrawMolecule(mol)
    compound_image.FinishDrawing()
    compound_image = compound_image.GetDrawingText()
    with open("../tmpimages/molecule.svg", "w") as f:
        f.write(compound_image)
