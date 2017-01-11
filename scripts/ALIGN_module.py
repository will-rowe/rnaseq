#!/usr/bin/env python
__author__ = 'Will Rowe'
__mail__ = "will.rowe@liverpool.ac.uk"

import os
import sys
import time
import subprocess
import glob
import GLOBALS # contains all the customisable variables that are used by the modules + main pipeline script
import logging
logger = logging.getLogger(__name__)

####
# Information
####
""" Alignment module

this module maps the supplied fastq files to the bowtie2 reference
"""


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
        logging.error('**** PIPELINE ERROR ****')
        logger.error('failed to create alignment directory')
        logger.error(e)
        sys.exit(1)

    # bowtie2
    logging.info('now running bowtie2 . . .')
    bt2_ref = '{}/temp_dir/reference_fasta.fa' .format(args.results_dir)
    bt2_cmd = 'bowtie2 -x {} -q {{}} --very-sensitive-local -p {} 2> {}/{{/.}}.bt2log | samtools view -@ {} -q {} -bS - | samtools sort -@ {} -m 1G - -o {}/{{/.}}.sorted.bam && samtools index {}/{{/.}}.sorted.bam' .format(bt2_ref, str(args.threads), alignment_dir, str(args.threads), GLOBALS.MAPQ_score, str(args.threads), alignment_dir, alignment_dir)
    parallel_bt2_cmd = 'printf \'{}\' | parallel -S {} --env PATH --workdir $PWD -j {} --delay 1.0 \'{}\'' .format('\\n'.join(sample_list), GLOBALS.SSH_list, GLOBALS.parallel_jobs, bt2_cmd)

    # run subprocess
    processes = []
    with open('{}/bt2.errorlog' .format(alignment_dir), 'a') as bt2_log:
        p1 = subprocess.Popen(parallel_bt2_cmd, shell=True, stdout=bt2_log)
        processes.append(p1)

    # wait for bt2 subprocess to complete:
    exit_codes = [p.wait() for p in processes]

    # count mapped reads (after MAPQ filter), remove samples from further analysis if no reads mapped
    bam_files = []
    for bam_file in glob.glob('{}/*.sorted.bam' .format(alignment_dir)):
        # note that this cmd only works for single end data
        samtools_cmd = 'samtools view -@ {} -F 0x904 -c {}' .format(str(args.threads), bam_file)
        try:
            p = subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE)
            init_mapped_reads, err = p.communicate()
        except Exception, e:
            logging.error('**** PIPELINE ERROR ****')
            logger.error('failed to generate bam files')
            logger.error(e)
            sys.exit(1)
        # flag warning if no reads align and remove sample
        if int(init_mapped_reads) == 0:
            logger.warning('No reads mapped for sample {}' .format(bam_file))
        # flag warning if fewer than 5 million reads align (arbitrary value)
        elif int(init_mapped_reads) <= 5000000:
            logger.warning('Low number ({}) of reads mapped for sample {}' .format(init_mapped_reads.rstrip("\n"), bam_file))
            bam_files.append(bam_file)
        else:
            bam_files.append(bam_file)

    # run multiqc once bowtie2 completes
    logging.info('generating multiqc report . . .')
    multiqc_cmd = 'multiqc --outdir {} --filename 00_multiqc_report {}' .format(args.results_dir, args.results_dir)
    processes = []
    with open('{}/00_multiqc.errorlog' .format(alignment_dir), 'a') as alignment_log:
        p1 = subprocess.Popen(multiqc_cmd, shell=True, stdout=alignment_log, stderr=alignment_log)
        processes.append(p1)

    # wait for multiqc subprocess to complete:
    exit_codes = [p.wait() for p in processes]

    return bam_files
