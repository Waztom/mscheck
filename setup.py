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
    install_requires=[
        "pyopenms==2.7.0",
        "scipy==1.9.3",
        "matplotlib==3.6.1",
        "svgutils==0.3.4",
    ],
    packages=find_packages(),
    url="https://github.com/Waztom/mscheck",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
