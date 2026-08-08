#!/usr/bin/env python
import io
import os
import re

from setuptools import setup, find_packages

version_match = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    io.open("mscheck/__init__.py", encoding="utf_8_sig").read(),
)
if version_match is None:
    raise ValueError("Version could not be determined")
__version__ = version_match.group(1)

if os.path.exists("README.md"):
    long_description = open("README.md").read()
else:
    long_description = """Auto MS mass checker"""

setup(
    name="mscheck",
    version=__version__,
    author="Warren Thompson",
    author_email="waztom@gmail.com",
    py_modules=["mscheck"],
    description="Auto MS mass checker",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.10",
    install_requires=[
        "rdkit>=2022.9",
        "matplotlib>=3.3",
        "scipy>=1.7",
        "numpy>=1.21",
        "svgutils>=0.3",
        "plotly>=5.0",
        "pandas>=1.3",
        "seaborn>=0.11",
        "PyYAML>=5.4",
        "numpyencoder>=0.3",
        "rainbow-api>=1.0",
        "psutil>=5.9",
    ],
    packages=find_packages(),
    url="https://github.com/xchem/mscheck",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
