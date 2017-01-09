#!/usr/bin/env python2
__author__ = 'Will Rowe'
__mail__ = "will.rowe@liverpool.ac.uk"

# Install a project in editable mode
# pip install -e . --user

from setuptools import setup

setup(
    name = "rnaseq",
    packages = ["rnaseq"],
    package_dir = {'rnaseq': 'scripts'},
    scripts = ['scripts/dupRadar.sh'
    ],
    entry_points = {
        "console_scripts": ['rnaseq = rnaseq.rnaseq:main']
        },
    version = '0.1.0',
    description = 'RNA-seq analysis pipeline for the Hinton Lab, Liverpool, UK',
    long_description = open('README.md').read(),
    author = 'Will Rowe',
    author_email = 'will.rowe@liverpool.ac.uk',
    url = 'http://will-rowe.github.io/',
    license = 'license.txt',
    install_requires = [
        "numpy == 1.7.0",
        "matplotlib == 1.5.3",
        "Jinja2 == 2.7.3",
        "MarkupSafe == 0.23",
        "multiqc == 0.8",
    ],
)
