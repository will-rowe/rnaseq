#!/usr/bin/env python
__author__ = 'Will Rowe'
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
#SSH_list = 'ada04,ada05,ada06,ada07,ada08,ada09,ada10,ada11'
parallel_jobs = '1'

# Settings for QC
trimmomatic_min_length = '40'
adapter_file = '/pub46/willr/000_HOME/0003_SOFTWARE/rnaseq/demo-files/TruSeq3-SE.fa'
kraken_db = '/pub46/willr/000_HOME/0002_REF_DATA/0004_DBs/minikraken_20141208'
#kraken_db = '/Users/willrowe/Documents/kraken/minikraken_20141208'
kraken_search_term = 'Salmonella\n' # include \n to only match Genus in Kraken report file

# Settings for alignment
MAPQ_score = '10'
