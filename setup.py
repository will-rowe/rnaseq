#!/usr/bin/env python2
__author__ = 'Will Rowe'
__mail__ = "will.rowe@liverpool.ac.uk"

# Install a project in editable mode
# pip install -e . --user

from setuptools import setup

setup(
    name = "rnaseq",
    version = '0.1.0',
    packages = ["rnaseq"],
    author = 'Will Rowe',
    author_email = 'will.rowe@liverpool.ac.uk',
    url = 'http://will-rowe.github.io/',
    description = 'RNA-seq analysis pipeline for the Hinton Lab, Liverpool, UK',
    long_description = open('README.md').read(),
    package_dir = {'rnaseq': 'scripts'},
    scripts = ['scripts/dupRadar.sh'
    ],
    entry_points = {
        "console_scripts": ['rnaseq = rnaseq.rnaseq:main']
        },
    license = 'license.txt',
    install_requires = [
        "numpy == 1.7.0",
        "matplotlib == 1.5.3",
        "Jinja2 == 2.7.3",
        "MarkupSafe == 0.23",
        # the commandline executable multiqc needs to be installed and in path
        #"multiqc == 0.8",
    ],
)
