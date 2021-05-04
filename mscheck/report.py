"""Genearte report function"""
from __future__ import annotations
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
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
def create_report_plot(
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
        fig, ax = plt.subplots(figsize=(12, 4))
        fig.suptitle(
            "MSCheck report: Mass not found for {}".format(compound_name), size=16
        )

        ax.plot(RT_values, TIC_values, color="k", zorder=-1)
        ax.set_xlabel("Retention time (min)")
        ax.set_ylabel("Total ion count (TIC)")
        ax.set_xlim([-0.9, ax.get_xlim()[1]])
        ax.set_ylim([0, ax.get_ylim()[1]])

        fig.savefig("../tmpimages/plot.svg", transparent=True)

        create_molecule_svg(mol)

        Figure(
            "29cm",
            "40cm",
            SVG("../tmpimages/molecule.svg").scale(0.004).move(3, 3),
            SVG("../tmpimages/plot.svg").scale(0.03),
        ).save("../reports/{}-report.svg".format(compound_name))

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
        ax[0].set_xlim([-0.9, ax[0].get_xlim()[1]])
        ax[0].set_ylim([0, ax[0].get_ylim()[1]])

        subplot = 1

        colors = iter(cm.rainbow(np.linspace(0, 1, no_plots * 2)))
        for ion_found, RT_values, TIC_values, mz_data, mz_strongest in zip(
            match_data["ions"],
            match_data["RT"],
            match_data["TIC"],
            match_data["mz_data"],
            match_data["mz_strongest"],
        ):

            ion_name = ion_found[0].strip("[]")
            mz_masses_max, mz_intensities_max, max_index = mz_strongest

            color_matches = next(colors)
            ax[0].scatter(
                [RT for i, RT in enumerate(RT_values) if i != max_index],
                [TIC for i, TIC in enumerate(TIC_values) if i != max_index],
                color=color_matches,
                s=45.0,
                linewidth=3,
                marker="x",
                zorder=1,
                label="{} ion matches".format(ion_name),
            )
            ax[0].legend(loc="upper right")

            RT_max = RT_values[max_index]
            TIC_max = TIC_values[max_index]
            color_max = next(colors)
            ax[0].scatter(
                RT_max,
                TIC_max,
                color=color_max,
                s=45.0,
                linewidth=3,
                marker="o",
                zorder=1,
                label="Strongest mz match".format(ion_name),
            )
            ax[0].legend(loc="upper right")

            markerline, stemline, baseline = ax[subplot].stem(
                mz_masses_max, mz_intensities_max
            )
            plt.setp(markerline, "markerfacecolor", color_max, "markersize", 10)
            plt.setp(stemline, "color", "k")
            plt.setp(baseline, "color", "k")

            ax[subplot].set_title(
                "Stongest mz pattern matching M + {} ({}) at RT: {} (min) ".format(
                    ion_name, ion_found[1], np.round(RT_max, 1)
                )
            )
            ax[subplot].set_xlabel("m/z (Da)")
            ax[subplot].set_ylabel("Relative intensity")
            ax[subplot].set_xlim(
                [ax[subplot].get_xlim()[0], ax[subplot].get_xlim()[1] + 10]
            )

            for i, j in zip(mz_masses_max, mz_intensities_max):
                annotation = ax[subplot].annotate(
                    str(i),
                    xy=(i, j),
                    textcoords="offset points",
                    xytext=(5.5, -3.5),
                    ha="left",
                )

            subplot += 1

        fig.savefig(
            "../tmpimages/plot.svg",
            transparent=True,
        )

        create_molecule_svg(mol)

        Figure(
            "29cm",
            "40cm",
            SVG("../tmpimages/molecule.svg").scale(0.004).move(3, 3),
            SVG("../tmpimages/plot.svg").scale(0.03),
        ).save("../reports/{}-report.svg".format(compound_name))
