from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import ntpath
import os
import pandas as pd
from .analyse import AnalyseSpectrum


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


def batch_analyse(
    csv_input_path: str,
    csv_output_path: str,
    analysis_types: list,
    data_dir: str,
    report_dir: str,
    modes: list,
) -> None:
    """
    Bulk analyse using a csv file containing target compounds and other metadata
    Args:
        csv_input_path (str): path to csv file containing batch data to analyse
        csv_output_path (str): path to save output csv file
        analysis_types (list): list of analysis types to perform eg. product, intermediate, reactant
        data_dir (str): path to the directory containing mzML files
        report_dir (str): path to the directory to save reports to
        modes (list): list of ionisation modes to analyse
    """

    batch_data = pd.read_csv(csv_input_path)

    for batch_index, batch_row in batch_data.iterrows():
        mzML_filename = batch_row["mzML-filename"]
        mzML_filepath = os.path.join(data_dir, mzML_filename, ".mzML")

        for analysis_type in analysis_types:
            no_analysis_type = batch_row["no-{}s".format(analysis_type)]
            analysis_type_ions_to_add = batch_row[
                "{}-ions-to-add".format(analysis_type)
            ].split(",")
            analysis_type_ions_to_sub = batch_row[
                "{}-ions-to-sub".format(analysis_type)
            ].split(",")
            analysis_type_match_tolerance = batch_row[
                "{}-match-tolerance".format(analysis_type)
            ]
            for no_analysis_type in range(no_analysis_type):
                analysis_type_smiles = batch_row[
                    "{}-{}".format(analysis_type, no_analysis_type + 1)
                ]
                for mode in modes:
                    analysis_type_analysis = AnalyseSpectrum(
                        mzMLfilepath=mzML_filepath, mode=mode
                    )
                    analysis_type_analysis.analyse(
                        compoundsmiles=analysis_type_smiles,
                        ionstoadd=analysis_type_ions_to_add,
                        ionstosub=analysis_type_ions_to_sub,
                        tolerance=analysis_type_match_tolerance,
                    )
                    analysis_type_analysis.create_report(
                        folder=os.path.join(report_dir, mode)
                    )
                    batch_row[
                        "{}-{}-max-EIC-signal".format(
                            analysis_type, no_analysis_type + 1
                        )
                    ] = analysis_type_analysis.analysedata["max_EIC_signal"]
                    batch_row[
                        "{}-{}-max-mz-match".format(analysis_type, no_analysis_type + 1)
                    ] = analysis_type_analysis.analysedata["max_mz_match"]
                    batch_row[
                        "{}-{}-ions-matched".format(analysis_type, no_analysis_type + 1)
                    ] = analysis_type_analysis.analysedata["ions"]

    batch_data.to_csv(csv_output_path, index=False)
