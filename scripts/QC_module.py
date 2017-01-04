#!/usr/bin/env python
__author__ = 'Will Rowe'
__version__ = 1.0
__mantainer__ = "Will Rowe"
__mail__ = "will.rowe@liverpool.ac.uk"

import os
import sys
import time
from math import sqrt
import subprocess
import glob
import re
import GLOBALS # contains all the customisable variables that are used by the modules + main pipeline script
import logging
logger = logging.getLogger(__name__)

####
# Information
####
""" QC module

this module quality checks all of the samples (using FastQC, Kraken, Trimmomatic and Dupradar)
flags are raised if samples look suspect
QC report generated using multiqc
"""


####
# Functions
####
def run_SD(input_values, population=True):
    """ Function to calculate the standard deviation from the mean kraken percentage
    reads in a dictionary (key = file name, value = percentage), raises flag for any sample outside of 1 SD from mean
    this has been adapted from codeselfstudy.com
    """
    num_items = len(input_values)
    mean = sum(input_values.values()) / num_items
    differences = [x - mean for x in input_values.values()]
    sq_differences = [d ** 2 for d in differences]
    ssd = sum(sq_differences)

    if population is True:
        # POPULATION standard deviation
        variance = ssd / num_items
    else:
        # SAMPLE standard deviation
        variance = ssd / (num_items - 1)
    sd = sqrt(variance)

    # flag warnings for any sample with Kraken result for keyword (Salmonella) outside of 1 SD of mean
    top = mean + sd
    bottom = mean - sd
    for key in input_values:
        if input_values[key] > float(top) or input_values[key] < float(bottom):
            logger.warning(' * a sample has kraken results outside 1 SD of collective mean ({})' .format(key))


def run_initQC(args, sample_list):
    """ Function to run an initial fastqc check, contamination screen and quality trim on all samples
    fastqc and kraken are run first - kraken is used to check that all samples have similar levels of our taxa of interest
    trimmomatic is used to quality trim the reads and remove adapters
    once QC has run, multiqc generates a report for all samples
    """
    # set up initial QC
    QC_dir = '{}/QC' .format(args.results_dir)
    try:
        os.makedirs(QC_dir)
    except Exception, e:
        run_sendError(e)
    half_threads = int(args.threads) / 2

    # initial fastqc run
    initial_FASTQC_dir = '{}/initial_FASTQC_results' .format(QC_dir)
    os.makedirs(initial_FASTQC_dir)
    fastqc_cmd = 'fastqc --threads {} --quiet --outdir {} {{}}' .format(str(half_threads), initial_FASTQC_dir)
    parallel_fastq_cmd = 'printf \'{}\' | parallel -S {} --env PATH --workdir $PWD -j {} --delay 1.0 \'{}\'' .format('\\n'.join(sample_list), GLOBALS.SSH_list, GLOBALS.parallel_jobs, fastqc_cmd)

    # initial kraken screen
    initial_KRAKEN_dir = '{}/initial_KRAKEN_results' .format(QC_dir)
    os.makedirs(initial_KRAKEN_dir)
    kraken_cmd = 'kraken --threads {} --preload --fastq-input --gzip-compressed --db {} {{}} | kraken-report --db {} > {}/{{/.}}-krakenreport.txt' .format(str(half_threads), GLOBALS.kraken_db, GLOBALS.kraken_db, initial_KRAKEN_dir)
    parallel_kraken_cmd = 'printf \'{}\' | parallel -S {} --env PATH --workdir $PWD -j {} --delay 1.0 \'{}\'' .format('\\n'.join(sample_list), GLOBALS.SSH_list, GLOBALS.parallel_jobs, kraken_cmd)

    # run subprocesses
    processes = []
    with open('{}/fastqc.errorlog' .format(initial_FASTQC_dir), 'a') as fastqc_log:
        p1 = subprocess.Popen(parallel_fastq_cmd, shell=True, stdout=fastqc_log, stderr=fastqc_log)
        processes.append(p1)
    with open('{}/kraken.errorlog' .format(initial_KRAKEN_dir), 'a') as kraken_log:
        p2 = subprocess.Popen(parallel_kraken_cmd, shell=True, stdout=kraken_log, stderr=kraken_log)
        processes.append(p2)

    # wait for fastqc and kraken subprocesses to complete:
    exit_codes = [p.wait() for p in processes]

    # use krakenreport to find number of salmonella sequences (this section of code currently discards secondary hits to the keyword (Salmonella\n))
    logging.info('collecting QC results . . .')
    kraken_percentages = {}
    for kraken_report in glob.glob('{}/*-krakenreport.txt' .format(initial_KRAKEN_dir)):
        with open(kraken_report, 'r') as fh:
            krakenreport_lines = fh.readlines()
            kraken_search = [line for line in krakenreport_lines if re.search(GLOBALS.kraken_search_term, line)]
            if kraken_search:
                kraken_percentages[kraken_report] = float(kraken_search[0].split('\t')[0])
            else:
                logging.error(' * no taxa hits to key word ({}) in sample {}' .format(GLOBALS.kraken_search_term.rstrip('\n'), kraken_report))

    # calculate mean and standard deviation for kraken percentages, then flag any samples outside 1 SD of the mean
    if len(kraken_percentages) < 1:
        logging.error(' * kraken found no hits for the keyword {} in any sample' .format(GLOBALS.kraken_search_term.rstrip('\n')))
    else:
        run_SD(kraken_percentages)

    # trimmomatic
    trimmomatic_dir = '{}/trimmomatic_data' .format(QC_dir)
    os.makedirs(trimmomatic_dir)
    trim_cmd = 'trimmomatic SE -threads {} {{}} {}/{{/.}}.trimmed.fq SLIDINGWINDOW:4:20 MINLEN:{} &> {}/{{/.}}-trimmomatic.log' .format(str(args.threads), trimmomatic_dir, GLOBALS.trimmomatic_min_length, trimmomatic_dir)
    parallel_trim_cmd = 'printf \'{}\' | parallel -S {} --env PATH --workdir $PWD -j {} --delay 1.0 \'{}\'' .format('\\n'.join(sample_list), GLOBALS.SSH_list, GLOBALS.parallel_jobs, trim_cmd)

    # run subprocess
    processes = []
    with open('{}/trimmomatic.errorlog' .format(trimmomatic_dir), 'a') as trim_log:
        p1 = subprocess.Popen(parallel_trim_cmd, shell=True, stdout=trim_log, stderr=trim_log)
        processes.append(p1)

    # wait for trimmomatic subprocesses to complete:
    exit_codes = [p.wait() for p in processes]

    # run multiqc once QC completes
    logging.info('generating multiqc report . . .')
    multiqc_cmd = 'multiqc --outdir {} --filename 00_multiqc_report {}' .format(QC_dir, QC_dir)
    processes = []
    with open('{}/00_multiqc.errorlog' .format(QC_dir), 'a') as qc_log:
        p1 = subprocess.Popen(multiqc_cmd, shell=True, stdout=qc_log, stderr=qc_log)
        processes.append(p1)

    # wait for trimmomatic subprocesses to complete:
    exit_codes = [p.wait() for p in processes]

    # return a list of trimmed files:
    trimmed_files = []
    for trimmed_file in glob.glob('{}/*trimmed.fq' .format(trimmomatic_dir)):
        trimmed_files.append(trimmed_file)
    return trimmed_files


def run_dupRadar(GFF_file, args):
    """ Function to run dupRadar.R script

    this function will be added to the pipeline at a later date
    currently there are no checks for the dupRadar library
    """

    dupradar_cmd = 'dupRadar.R infile.bam {} stranded=yes paired=no outdir=./dupRadar threads={}' .format(GFF_file.gff_filename, str(args.threads))
    print dupradar_cmd
