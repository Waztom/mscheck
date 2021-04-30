import matplotlib.pyplot as plt
from analyse import AnalyseSpectrum


class ReportSpectrum(AnalyseSpectrum):
    """
    Creates a report of an analysed spectrum
    """


def get_MS_plot(self):
    """
    Creates a plot of the total ion count spectrum
    """
    plt.scatter(self.MSdata["RT"], self.MSdata["TIC"], s=0.5)
    plt.xlabel("Retention time (s)")
    plt.ylabel("TIC")
    plt.title(self._filepath)
    return plt.show()
