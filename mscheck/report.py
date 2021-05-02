from __future__ import annotations
import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, AnnotationBbox
import matplotlib.pylab as pylab
from svgutils.compose import *
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from pathlib import Path
from functools import wraps
import shutil

params = {
    "font.weight": "bold",
    "legend.fontsize": "x-large",
    "figure.figsize": (15, 5.85),
    "axes.labelsize": "large",
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "large",
    "ytick.labelsize": "large",
}
pylab.rcParams.update(params)


def create_report_dir(func):
    def wrapper(*args, **kwargs):
        Path("../reports").mkdir(parents=True, exist_ok=True)
        Path("../tmpimages").mkdir(parents=True, exist_ok=True)
        func(*args, **kwargs)
        shutil.rmtree("../tmpimages")

    return wrapper


def create_molecule_svg(mol: rdkitmol):
    """
    Creates svg image of rdkit mol and saves to file
    """
    compound_image = rdMolDraw2D.MolDraw2DSVG(824, 556)
    compound_image.DrawMolecule(mol)
    compound_image.FinishDrawing()
    compound_image = compound_image.GetDrawingText()
    with open("../tmpimages/molecule.svg", "w") as f:
        f.write(compound_image)


def get_strongest_mz_pattern(mz_values: list) -> largest_mz_pattern:
    """
    Finds largest mz signal from list of mz patterns
    """
    mz_intensities = [mz_value[1] for mz_value in mz_values]
    mz_masses = [mz_value[0] for mz_value in mz_values]

    for index, intensity_values in enumerate(mz_intensities):
        if index == 0:
            max_value = np.amax(intensity_values)
            max_index = index
        else:
            max_value_test = np.amax(intensity_values)
            if max_value_test > max_value:
                max_value = max_value_test
                max_index = index
            else:
                pass
    mz_masses_max = mz_masses[max_index]
    mz_intensities_max = mz_intensities[max_index]
    return mz_masses_max, mz_intensities_max, max_index


@create_report_dir
def create_plot_report(
    RT_values: list,
    TIC_values: list,
    plot_title: str,
    no_plots: int,
    mol: rdkitmol,
    match_data: dict = None,
) -> plot:
    if no_plots == 0:
        fig, ax = plt.subplots()
        ax.plot(RT_values, TIC_values, color="k", zorder=-1)
        ax.set_xlabel("Retention time (min)")
        ax.set_ylabel("Total ion count (TIC)")
        ax.set_title(plot_title)
        fig.savefig("../tmpimages/plot.svg", transparent=True)

        create_molecule_svg(mol)

        Figure(
            "31.0cm",
            "22.0cm",
            SVG("../tmpimages/molecule.svg").scale(0.012).move(16.5, 1.5),
            SVG("../tmpimages/plot.svg").scale(0.05),
        ).save("../reports/{}report.svg".format(plot_title))
    else:
        fig, ax = plt.subplots(no_plots + 1)
        ax[0].plot(RT_values, TIC_values, color="k", zorder=-1)
        ax[0].set_xlabel("Retention time (min)")
        ax[0].set_ylabel("Total ion count (TIC)")
        ax[0].set_title(plot_title)

        subplot = 1

        for ion_found, RT_values, TIC_values, mz_values in zip(
            match_data["ions_matched"],
            match_data["RT_matched"],
            match_data["TIC_matched"],
            match_data["mz_patterns"],
        ):

            ax[0].scatter(RT_values, TIC_values, s=25.0, marker="x", color="r", zorder=1)
            TIC_max_index = np.argmax(TIC_values)

            mz_masses_max, mz_intensities_max, max_index = get_strongest_mz_pattern(
                mz_values
            )

            ax[subplot].stem(mz_masses_max, mz_intensities_max)
            ax[subplot].set_title(
                "Stongest mz pattern matching M + {} = {}".format(
                    ion_found[0], ion_found[1]
                )
            )

            for i, j in zip(mz_masses_max, mz_intensities_max):
                ax[subplot].annotate(
                    str(i),
                    xy=(i, j),
                    textcoords="offset points",
                    xytext=(0, 1),
                    ha="left",
                )

            xy = (RT_values[max_index], TIC_values[max_index])

            offsetbox = TextArea(
                "Ion: {} \nM+: {}\nRT: {}".format(
                    ion_found[0], ion_found[1], np.round(RT_values[TIC_max_index], 2)
                )
            )

            ab = AnnotationBbox(
                offsetbox,
                xy,
                xybox=(20, 50),
                xycoords="data",
                boxcoords=("offset points"),
                box_alignment=(0.3, 0.4),
                arrowprops=dict(arrowstyle="->"),
            )
            ax[0].add_artist(ab)

            subplot += 1

        plt.subplots_adjust(
            top=2.0, bottom=0.01, left=0.10, right=0.95, hspace=0.3, wspace=0.5
        )

        fig.savefig("../tmpimages/plot.svg", transparent=True)

        create_molecule_svg(mol)

        Figure(
            "100.0cm",
            "100.0cm",
            SVG("../tmpimages/molecule.svg").scale(0.012).move(16.5, 1.5),
            SVG("../tmpimages/plot.svg").scale(0.3),
        ).save("../reports/{}-report.svg".format(plot_title))
