# <a name="MScheck for locating target compound masses in mass spectra"></a>**MScheck for locating target compound masses in mass spectra**

[//]: # "Badges"

[![PyPI version](https://badge.fury.io/py/mscheck.svg)](https://badge.fury.io/py/mscheck)

MScheck is a python package that hunts for a target compound mass + given ion mass (eg. H+, Na+).
MScheck was created to assist with the automated mass spectrum analysis of target compounds synthesised using
a high throughput approach.

# <a name="MScheck installation"></a>**MScheck installation**

MScheck relies on rdkit for generating molecule SVG images and for calculating molecular weights<br>
Installing rdkit using conda works best followed by a pip install of MScheck<br>

1. Create a conda environment with rdkit

   > `conda create -c conda-forge -n MScheck rdkit`

2. Activate the MScheck conda environment created

   > `conda activate MScheck`

3. Pip install MScheck
   > `pip install mscheck`

# <a name="MScheck use"></a>**MScheck use**

MScheck has been tested on using Agilent LCMS files (.D) as the starting file format. One challenge is to convert vendor file formats into
a format with the binary decoded.

1. Convert .D Agilent folder to a .d MassHunter format using Agilent's ChemStation to MassHunter Translator (B.04.00)
2. Convert .d format into .mzML format using [ProteoWizard's](http://proteowizard.sourceforge.net/) MSConvert tool
3. Finally - we have an file format that we can use!
4. In your favourite IDE or Jupyter notebook - as basic example of using MSCheck is prvided below:<br>

```
from mscheck.analyse import AnalyseSpectrum

# Create MS scptrum object and find peaks
test = AnalyseSpectrum("<path to .mzML file>", mode="Positive")

# Analyse test spectrum
test.analyse(compoundsmiles=target_compound,
             ionstoadd=["[H]", "[Na]", "[K]", "[NH4+]"],
             tolerance=1)test.analyse

# Create a SVG report - if you do not give a compound_name
# the ending leaf of the file name will be used
test.create_report(compound_name="Test")
```

 <br>

Example of report output:<br>

<p align="center">
<img src="images/report.svg" width="600px">
</p>

The .SVG report will be in a folder called Reports
