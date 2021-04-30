#!/usr/bin/env python
import io
import os
import re

from setuptools import setup, find_packages
version_match = re.search(
    r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
    io.open('pycule/__init__.py', encoding='utf_8_sig').read()
)
if version_match is None:
    raise ValueError('Version could not be determined')
__version__ = version_match.group(1)

if os.path.exists('README.md'):
    long_description = open('README.md').read()
else:
    long_description = '''Python wrapper for the Ultimate MCule API'''

setup(
    name='PyCule',
    version=__version__,
    author='Warren Thompson',
    author_email='waztom@gmail.com',
    py_modules=['PyCule'],
    description='Python wrapper for the MCule Ultimate API',
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    install_requires=[
        'requests==2.23.0',
    ],
    packages=find_packages(),
    url='https://github.com/Waztom/pycule',
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)