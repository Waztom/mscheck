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

from utils import create_molecule_svg

params = {
    "font.weight": "bold",
    "legend.fontsize": "small",
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
        target_folder = kwargs["folder"]
        Path("../{}".format(target_folder)).mkdir(parents=True, exist_ok=True)
        Path("../tmpimages").mkdir(parents=True, exist_ok=True)
        func(*args, **kwargs)
        shutil.rmtree("../tmpimages")

    return wrapper


@create_report_dir
def create_report_plot(
    msmode: str,
    RT_values: list,
    TIC_values: list,
    compound_name: str,
    no_plots: int,
    mol: rdkitmol,
    folder: str,
    analysedata: dict = None,
) -> plot:
    """
    Creates report
    """
    if msmode == "Positive":
        ion_mode = "+"
    if msmode == "Negative":
        ion_mode = "-"
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
        ).save("../{}/{}-report.svg".format(folder, compound_name))
        plt.close("all")

    else:
        fig, ax = plt.subplots((no_plots * 2) + 1)

        fig.tight_layout(pad=4.2)

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
        for EIC_data, max_mz_match, ion_found, RT_match_values, TIC_match_values, mz_strongest in zip(
            analysedata["EIC_data"],
            analysedata["max_mz_match"],
            analysedata["ions"],
            analysedata["RT"],
            analysedata["TIC"],
            analysedata["mz_strongest"],
        ):
            if not max_mz_match:
                max_match_label = "Not dominant"
            else:
                max_match_label = "Dominant"

            ion_name = ion_found[0].strip("[]")
            mz_masses_max, mz_intensities_max, max_index = mz_strongest

            color_matches = next(colors)              

            ax[0].scatter(
                [RT for i, RT in enumerate(RT_match_values) if i != max_index],
                [TIC for i, TIC in enumerate(TIC_match_values) if i != max_index],
                color=color_matches,
                s=45.0,
                linewidth=3,
                marker="x",
                zorder=1,
                label="{} mz match for M{}{}".format(max_match_label, ion_mode, ion_name),
            )
            ax[0].legend(loc="upper right")

            RT_max = RT_match_values[max_index]
            TIC_max = TIC_match_values[max_index]
            color_max = next(colors)
            ax[0].scatter(
                RT_max,
                TIC_max,
                color=color_max,
                s=45.0,
                linewidth=3,
                marker="o",
                zorder=1,
                label="Strongest mz match pattern for M{}{}".format(ion_mode,ion_name),
            )
            ax[0].legend(loc="upper right")

            markerline, stemline, baseline = ax[subplot].stem(
                mz_masses_max, mz_intensities_max
            )
            plt.setp(markerline, "markerfacecolor", color_max, "markersize", 10)
            plt.setp(stemline, "color", "k")
            plt.setp(baseline, "color", "k")

            ax[subplot].set_title(
                "Stongest mz pattern matching M{}{} ({}) at RT: {} (min) ".format(
                    ion_mode,ion_name, ion_found[1], np.round(RT_max,1)
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

            ax[subplot].set_title("EIC for M{}{}".format(ion_mode,ion_name))
            ax[subplot].set_xlabel("Retention time (min)")
            ax[subplot].set_ylabel("Extracted ion intensity")
            ax[subplot].set_xlim([-0.9, ax[0].get_xlim()[1]])
    
            ax[subplot].plot(
                [data[0] for data in EIC_data],
                [data[1] for data in EIC_data],
                color=color_max,
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
            SVG("../tmpimages/molecule.svg").scale(0.0035).move(3.4, 2.0),
            SVG("../tmpimages/plot.svg").scale(0.03),
        ).save("../{}/{}-report.svg".format(folder, compound_name))
        plt.close("all")
