# <a name="MScheck for locating target compound masses in mass spectra"></a>**MScheck for locating target compound masses in mass spectra**

[//]: # "Badges"

[![PyPI version](https://badge.fury.io/py/mscheck.svg)](https://badge.fury.io/py/mscheck)

[![build test](https://github.com/waztom/mscheck/actions/workflows/build-test.yml/badge.svg)](https://github.com/waztom/mscheck/actions/workflows/build-test.yml)

MScheck is a python package that hunts for a target compound mass + given ion mass (eg. H+, Na+).
MScheck was created to assist with the automated mass spectrum analysis of target compounds synthesised using
a high throughput approach.

MSCheck supports two input formats:

| Format | Extension | MS | UV / DAD | ELSD |
|---|---|---|---|---|
| Open mzML | `.mzML` | ✅ | ✅ (single-wavelength + full 2D DAD when present) | ❌ |
| Agilent MassHunter | `.D` / `.dx` | ✅ | ✅ | ✅ |
| Waters MassLynx | `.raw` | ✅ | ✅ | ✅ |

mzML is the preferred open interchange format. Vendor directories (`.D`, `.raw`) are read directly via [rainbow-api](https://github.com/evanyeyeye/rainbow), which unlocks ELSD data that is not carried by mzML conversions.

The `AnalyseMS` class uses Scipy's signal peak finding algorithms (`find_peaks` and `peak_widths`) to find peaks and calculate the FWHM. Mass spectrum data points are analysed in the area of the peak above the FWHM height by searching for the sum of the parent mass of the target molecule and ion. Different ions can be included in the search — see the example below.

# <a name="MScheck installation"></a>**MScheck installation**

**Note:** `rdkit` installs cleanly from PyPI (`pip install rdkit`) on Python ≥ 3.10.
A virtual environment is strongly recommended to avoid dependency conflicts.

### Option A — Python venv

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install mscheck
```

### Option B — conda environment

```bash
conda create -n mscheck python=3.11
conda activate mscheck
pip install mscheck
```

# <a name="Input files"></a>**Input files**

MSCheck accepts:

- **`.mzML` files** — the open PSI standard; produced by [ProteoWizard msConvert](http://proteowizard.sourceforge.net/) from any vendor format. This is the recommended path for sharing and archiving data.
- **Agilent `.D` / `.dx` directories** — read directly without conversion. Provides full UV 2D DAD and ELSD data.
- **Waters `.raw` directories** — read directly without conversion. Provides UV and ELSD data.

To convert Agilent `.D` files to mzML (optional):

1. Convert `.D` Agilent folder to `.d` MassHunter format using Agilent's ChemStation to MassHunter Translator (B.04.00)
2. Convert `.d` format to `.mzML` using [ProteoWizard msConvert](http://proteowizard.sourceforge.net/)

In your favourite IDE or Jupyter notebook — a basic example of using MSCheck is provided below:<br>

```python
from mscheck.analyseMS import AnalyseMS

# Accepts .mzML files or Agilent .D / Waters .raw directories
test = AnalyseMS("<path to .mzML file or .D directory>", mode="Positive")

# Set SMILES of target to search for
target_SMILES = "CCOC(=O)N1CCN(C(=O)N2CCN(C(=O)c3ccco3)CC2)CC1"

# Analyse test spectrum searching for target SMILES
test.analyse(compoundsmiles=target_SMILES,
             ionstoadd=["[H]", "[Na]", "[K]", "[NH4+]"],
             tolerance=1)

# Create a .svg report
test.create_report(compound_name="Test")
```

 <br>

The .svg report will be in a folder called Reports

For a runnable, end-to-end walkthrough (including bulk analysis with `BulkAnalyser`), see [notebooks/mscheck_tutorial.ipynb](notebooks/mscheck_tutorial.ipynb).

Example of report output:<br>

<p align="center">
<img src="images/report.svg" width="600px">
</p>
