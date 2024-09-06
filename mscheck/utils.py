from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import ntpath
import os


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


def sort_dir_files(data_dir: str) -> list:
    """
    Sorts files by their numeric value
    """
    filelist = []

    for dir_, _, files in sorted(os.walk(data_dir)):
        #sorted(files)
        for file_name in sorted(files):
            rel_dir = os.path.relpath(dir_, data_dir)
            rel_file = os.path.join(rel_dir, file_name)
            if rel_file.endswith(".mzML"):
                filelist.append(rel_file)
    sorted_filelist = os_sorted(filelist)

    return sorted_filelist

def bulk_analyse(csv_input_path: str, data_dir: str):
        """
        Bulk analyse using a csv file containing target compounds and other metadata
        Args:
            csv_input_path (str): path to csv file containing target compounds and metadata
            data_dir (str): directory containing mzML files
        """

        target_data = pd.read_csv(csv_input_path)
        target_data["product-SMILES"] = target_data["product-SMILES"].astype(str)
        target_data["product-ions-to-add"] = target_data["product-ions-to-add"].astype(str)
        target_data["product-ions-to-sub"] = target_data["product-ions-to-sub"].astype(str)
        target_data["product-match-tolerance"] = target_data["product-match-tolerance"].astype(int)



        for index, row in target_data.iterrows():
            self.analyse(
                compoundsmiles=row["SMILES"],
                ionstoadd=row["Ionstoadd"].split(","),
                ionstosub=row["Ionstosub"].split(","),
                tolerance=row["Tolerance"],
            )
            self.create_report(folder=data_dir, compound_name=row["CompoundName"])
