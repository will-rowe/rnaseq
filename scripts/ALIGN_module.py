#!/usr/bin/env python
__author__ = 'Will Rowe'
__version__ = 1.0
__mantainer__ = "Will Rowe"
__mail__ = "will.rowe@liverpool.ac.uk"

import os
import sys
import time
import subprocess
import re
import GLOBALS # contains all the customisable variables that are used by the modules + main pipeline script
import logging
logger = logging.getLogger(__name__)

####
# Information
####
""" Alignment module

this module maps the supplied fastq files to the bowtie2 reference
it is a work in progress...
add in multiqc report
"""
# default settings



####
# Functions
####


def run_alignData(args, sample_list, GFF_file):
    """ Function to run Bowtie2
    returns a filtered bam file (MAPQ > 10)

    may not need GFF_file yet....
    """
    # set up alignment directory
    alignment_dir = '{}/ALIGNMENT' .format(args.results_dir)
    try:
        os.makedirs(alignment_dir)
    except Exception, e:
        run_sendError(e)


    # bowtie2
    bt2_ref = '{}/temp_dir/reference_fasta.fa' .format(args.results_dir)
    bt2_cmd = 'bowtie2 -x {} -q {{}} --very-sensitive-local -p {} | samtools view -q {} -bS - | samtools sort - -o {}/{{/.}}.sorted.bam && samtools index {}/{{/.}}.sorted.bam' .format(bt2_ref, str(args.threads), GLOBALS.MAPQ_score, alignment_dir, alignment_dir)
    parallel_bt2_cmd = 'printf \'{}\' | parallel -S {} --env PATH --workdir $PWD -j {} --delay 1.0 \'{}\'' .format('\\n'.join(sample_list), GLOBALS.SSH_list, GLOBALS.parallel_jobs, bt2_cmd)

    # run subprocess
    processes = []
    with open('{}/bt2.errorlog' .format(alignment_dir), 'a') as bt2_log:
        p1 = subprocess.Popen(parallel_bt2_cmd, shell=True, stdout=bt2_log, stderr=bt2_log)
        processes.append(p1)

    # wait for trimmomatic subprocesses to complete:
    exit_codes = [p.wait() for p in processes]
