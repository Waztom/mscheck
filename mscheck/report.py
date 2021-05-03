"""Genearte report function"""
from __future__ import annotations
import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, AnnotationBbox
import matplotlib.pylab as pylab
from svgutils.compose import *
import numpy as np
from pathlib import Path
import shutil

from .utils import create_molecule_svg

params = {
    "font.weight": "bold",
    "legend.fontsize": "x-large",
    "figure.figsize": (12, 12),
    "axes.labelsize": "large",
    "axes.titleweight": "bold",
    "axes.labelweight": "bold",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "large",
    "ytick.labelsize": "large",
    "figure.titleweight": "bold",
    "figure.subplot.bottom": 0.11,
    "figure.subplot.hspace": 0.2,
    "figure.subplot.left": 0.125,
    "figure.subplot.right": 0.9,
    "figure.subplot.top": 0.88,
    "figure.subplot.wspace": 0.2,
    "figure.constrained_layout.use": True,
}
pylab.rcParams.update(params)


def create_report_dir(func):
    def wrapper(*args, **kwargs):
        Path("../reports").mkdir(parents=True, exist_ok=True)
        Path("../tmpimages").mkdir(parents=True, exist_ok=True)
        func(*args, **kwargs)
        shutil.rmtree("../tmpimages")

    return wrapper


@create_report_dir
def create_plot_report(
    RT_values: list,
    TIC_values: list,
    compound_name: str,
    no_plots: int,
    mol: rdkitmol,
    match_data: dict = None,
) -> plot:
    """
    Creates report
    """

    if no_plots == 0:
        fig, ax = plt.subplots()
        fig.suptitle(
            "MSCheck report: Mass not found for {}".format(compound_name), size=16
        )

        ax.plot(RT_values, TIC_values, color="k", zorder=-1)
        ax.set_xlabel("Retention time (min)")
        ax.set_ylabel("Total ion count (TIC)")
        fig.savefig("../tmpimages/plot.svg", transparent=True)

        create_molecule_svg(mol)

        Figure(
            "31.0cm",
            "22.0cm",
            SVG("../tmpimages/molecule.svg").scale(0.012).move(16.5, 1.5),
            SVG("../tmpimages/plot.svg").scale(0.05),
        ).save("../reports/{}report.svg".format(compound_name))
    else:
        fig, ax = plt.subplots(no_plots + 1)

        fig.suptitle("MSCheck report: Mass found for {}".format(compound_name), size=16)

        fig.align_ylabels()

        ax[0].plot(
            RT_values,
            TIC_values,
            color="k",
            zorder=-1,
        )
        ax[0].set_xlabel("Retention time (min)")
        ax[0].set_ylabel("Total ion count (TIC)")
        ax[0].set_xlim([-0.7, ax[0].get_xlim()[1]])

        subplot = 1

        for ion_found, RT_values, TIC_values, mz_data, mz_strongest in zip(
            match_data["ions"],
            match_data["RT"],
            match_data["TIC"],
            match_data["mz_data"],
            match_data["mz_strongest"],
        ):

            ion_name = ion_found[0].strip("[]")
            mz_masses_max, mz_intensities_max, max_index = mz_strongest

            ax[0].scatter(
                RT_values,
                TIC_values,
                s=45.0,
                linewidth=3,
                marker="x",
                zorder=1,
                label="{} ion matches".format(ion_name),
            )
            ax[0].legend(loc="upper right")

            ax[subplot].stem(mz_masses_max, mz_intensities_max)
            ax[subplot].set_title(
                "Stongest mz pattern matching M + {} ({})".format(ion_name, ion_found[1])
            )
            ax[subplot].set_xlabel("m/z (Da)")
            ax[subplot].set_ylabel("Ion count")

            for i, j in zip(mz_masses_max, mz_intensities_max):
                annotation = ax[subplot].annotate(
                    str(i),
                    xy=(i, j),
                    textcoords="offset points",
                    xytext=(3.5, -3.5),
                    ha="left",
                )

            xy = (RT_values[max_index], TIC_values[max_index])

            offsetbox = TextArea(
                "Ion: {} \nM+: {}\nRT: {}".format(
                    ion_found[0], ion_found[1], np.round(RT_values[max_index], 2)
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

        fig.savefig(
            "../tmpimages/plot.svg",
            transparent=True,
        )

        create_molecule_svg(mol)

        test = Figure(
            "29cm",
            "40cm",
            SVG("../tmpimages/molecule.svg").scale(0.011).move(-1, 1.2),
            SVG("../tmpimages/plot.svg").scale(0.03),
        ).save("../reports/{}-report.svg".format(compound_name))
