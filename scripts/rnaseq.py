#!/usr/bin/env python2
__author__ = 'Will Rowe'
__version__ = 1.0
__mantainer__ = "Will Rowe"
__mail__ = "will.rowe@liverpool.ac.uk"

import os
import sys
import time
import argparse # used to parse the command line input
import logging # to create a pipeline log file
import gzip # used to handle gzipped files
import collections # used for OrderedDict
import getpass # for identifying username
import subprocess
import shutil # for removing temporary files
import GFF_module as GFF_module # required pipeline modules
import QC_module as QC_module # required pipeline modules
import GLOBALS as GLOBALS # contains all the customisable variables that are used by the modules + main script


####
# Information
####
""" The Hinton Lab RNA-seq pipeline



only works on single end illumina data
the pipeline checks for required programs but does not check for correct version numbers...
only checks file name extensions - could add more stringent fastq checks
add support for cluster and single compute node jobs
adapter removal
strain signature file


up to QC module - replacing fastqc with alternatives. Trying out multiqc report.
"""


####
# Classes
####
"""
class Analysis(object):

    Class to hold pipeline data and generate report


    def __init__(self, data_dir, plots_dir, samples):
        self.data_dir = data_dir
        self.plots_dir = plots_dir
        self.samples = samples

    def do_some_work(self):
        ...
        self.results = results

prj = Project("testprj")
prj.addSampleSheet("metadata/sample_annotation.csv")

analysis = Analysis("data", "results/plots", prj.samples)
analysis.do_some_work()
"""


class customParser(argparse.ArgumentParser):
    """ Class summary: a custom argument parser

    This class builds on the argparse parser to provide additional help by overiding the default error method
    """
    def __init__(self):
        description = 'The Hinton Lab RNA-seq pipeline',
        usage       = 'python rnaseq.py [OPTIONS] reference file_list'
        argparse.ArgumentParser.__init__(self, description, usage)

    def error(self, message):
        sys.stderr.write('error: {}\n\n' .format(message))
        self.print_help()
        sys.exit(2)


####
# Functions
####
def run_getArguments():
    # use our custom parser class that is built from the ArgumentParser (argparse library)
    parser = customParser()

    # set up positional arguments
    parser.add_argument('reference',
        help = 'reference to be used for read alignment / annotation --> must be in GFF3 format and include fasta (can be gzipped)')
    parser.add_argument('file_list',
        help = 'list of samples to be processed (newline separated, fastq files)')

    # set up optional arguments
    parser.add_argument('-c', '--checks', action = 'store_true', dest = 'checks',
        help = 'run setup checks (checks annotation is GFF3 formatted, bowtie2 indices are present etc.)', default = True)
    parser.add_argument('-o', '--outdir', default = './hinton-rnaseq-{}' .format(time.strftime('%H%M%S')), dest = 'results_dir',
        help = 'specify the directory to put the results files (default = ./hinton-rnaseq-xxxx)')
    parser.add_argument('-l', '--logfile', default = './log.txt', dest = 'log_file',
        help = 'specify the name of the log file')
    parser.add_argument('-t', '--threads', default = '10', dest = 'threads',
        help = 'specify the number of threads to use in the multithreading steps (default = 1)')
    parser.add_argument('-s', '--silent', action = 'store_true', dest = 'silent',
        help = 'silent behaviour', default = False)

    # parse the arguments and return them
    args = parser.parse_args()
    return args


def run_getLogger(args):
    # set up logging
    global logging
    logging.basicConfig(filename='{}/{}' .format(args.results_dir, args.log_file), level=logging.INFO, format='rnaseq.py: %(asctime)s\t%(levelname)s\t%(message)s', datefmt='%m/%d/%Y %I:%M:%S')
    return logging


def run_sendError(exception):
    logging.error('**** PIPELINE ERROR ****')
    logging.error(exception)
    sys.exit(1)


def run_pipelineSetup(args):
    # get the user
    user = getpass.getuser()

    # set up the output directory
    if not os.path.exists(args.results_dir):
        try:
            os.makedirs(args.results_dir)
        except Exception, e:
            run_sendError(e)

    # create the temp directory for pipeline
    temp_dir = '{}/temp_dir' .format(args.results_dir)
    try:
        os.makedirs(temp_dir)
    except Exception, e:
        run_sendError(e)

    # start the logger
    logger = run_getLogger(args)
    logging.info('RNA-seq pipeline started by {}' .format(user))
    logging.info('**** SETTING UP PIPELINE ****')

    # check required programs are installed
    logging.info('checking for required programs . . .')
    program_list = ['parallel', 'bowtie2', 'fastqc', 'kraken', 'kraken-report', 'trimmomatic', 'samtools', 'multiqc',]
    missing_programs = []
    for program in program_list:
        try:
            program_check = subprocess.check_output(['which', program])
        except subprocess.CalledProcessError:
            missing_programs.append(program)
            continue
        logging.info(' * {} found: {}' .format(program, program_check.rstrip()))
    if missing_programs:
        logging.error(' * can\'t find: {}' .format(','.join(missing_programs)))
        run_sendError('required programs not installed')

    # read the supplied annotation file
    logging.info('checking the input annotation file . . .')
    try:
        GFF_file = GFF_module.GFF_annotationFile(args.reference)
    except Exception, exception:
        run_sendError(exception)

    # check the annotation formatting, write a new gff file in the tmp_dir and then update our GFF_file object to point to our newly saved version
    checked_GFF_copy = '{}/reference_annotation.gff' .format(temp_dir)
    GFF_file.make_gffFile(checked_GFF_copy)
    GFF_file.gff_filename = checked_GFF_copy
    # extract the fasta
    logging.info('annotation file OK, extracting reference fasta . . .')
    reference_fasta = '{}/reference_fasta.fa' .format(temp_dir)
    with open(reference_fasta, 'w') as fasta:
        fasta.write(GFF_file.get_fastaSeqs())

    # setup the Bowtie2 index
    logging.info('creating bowtie2 index from reference fasta . . .')
    try:
        bowtie2_build = subprocess.check_output(['bowtie2-build', reference_fasta, reference_fasta], stderr=open('/dev/null', 'w'))
    except subprocess.CalledProcessError:
        run_sendError('can\'t run bowtie2-build command')

    # return the checked GFF file
    return GFF_file


def run_getInput(args):
    """ Simple input checker
    this funciton checks that the input file list exists, it then checks each file specified in the list before creating the sample run list for the pipeline
    """
    if not os.path.isfile(args.file_list):
        run_sendError('input file list doesn\'t appear to exist . . .')
    sample_file_list = [line.rstrip('\n') for line in open(args.file_list)]
    checked_samples = []
    checked_sample_file_list = '{}/sample_file_list.txt' .format(args.results_dir)
    logging.info('checking file list . . .')
    for i in sample_file_list:
        if not i.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
            logging.warning(' * skipping file with bad extension in file list: {}' .format(i))
        elif not os.path.isfile(i):
            logging.warning(' * pipeline can\'t find file: {}' .format(i))
        else:
            checked_samples.append(os.path.abspath(i))
    if len(checked_samples) == 0:
        run_sendError('no input files were passed to the pipeline . . .')
    else:
        logging.info(' * found {} fastq files present in input file list' .format(len(checked_samples)))
    return checked_samples


def main():
    # get the command line arguments
    args = run_getArguments()

    # set up the pipeline
    GFF_file = run_pipelineSetup(args)

    # get input files
    logging.info('**** COLLECTING INPUT FILES ****')
    sample_list = run_getInput(args)

    # QC the samples
    logging.info('**** RUNNING QC ****')
    trimmed_files = QC_module.run_initQC(args, sample_list)

    # align the samples to the reference genome





    # run dupRadar
    QC_module.run_dupRadar(GFF_file, args)






    # report pipeline finished and exit
    logging.info('**** FINISHED ****')
    logging.info("pipeline completed successfully")
    #shutil.rmtree(temp_dir)
    sys.exit(1)


####
# Main
####
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\nPipeline canceled by user!")
        sys.exit(0)
