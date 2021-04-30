from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from PIL import Image
from io import BytesIO
import matplotlib.pyplot as plt


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


def get_molecule_image(mol):
    """
    Creates a PIl image of a mol object
    Args:
        mol: mol object target compound
    """
    image = Draw.MolToImage(mol)
    figfile = BytesIO()
    image.save(figfile, bbox_inches="tight", format="png")
    image = Image.open(figfile)
    image = image.resize((500, 500))
    return image
