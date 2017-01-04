#!/usr/bin/env python
__author__ = 'Will Rowe'
__version__ = 1.0
__mantainer__ = "Will Rowe"
__mail__ = "will.rowe@liverpool.ac.uk"


####
# Information
####
""" Global variables

this file contains global variables for the Hinton Lab RNA-seq pipeline, such as the CGR nodes to use in the parallel steps etc.
"""


####
# Globals
####
# GNU Parallel settings
SSH_list = 'ada04,ada05,ada06,ada07,ada08,ada09,ada10,ada11'
parallel_jobs = '1'

# Settings for QC
trimmomatic_min_length = '40'
kraken_db = '/pub46/willr/000_HOME/0002_REF_DATA/0004_DBs/minikraken_20141208'
kraken_search_term = 'Salmonella\n' # include \n to only match Genus in Kraken report file
