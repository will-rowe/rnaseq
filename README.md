# RNA-seq pipeline
A pipeline for analysis of Salmonella RNA-seq data

-----------


## Contents

* [Overview](#overview)
* [To Do](#todo)
* [Installation](#installation)
* [Running](#running)
* [Output](#output)
* [Known Issues](#issues)


<a name="overview"/>

## Overview

This is a generic pipeline for processing RNA-seq data. It is intended to standardise the workflow used by the Hinton Lab (University of Liverpool, UK) to analyse stranded, single-end Illumina RNA-seq libraries. The software/workflow used is:

1. QC

  * **FastQC**: assess read quality/length, adapters etc.
  * **Kraken**: check for our typical taxonomic signature
  * **Trimmomatic**: remove adapters and quality trim reads
  * **DupRadar**: check for PCR artifacts

2. Alignment

  * **Bowtie2**: align reads to a reference genome
  * **Samtools**: alignment post-processing
  * **Bedtools**: are the usual transcripts covered? (QC)
  * **FastQC**: repeat initial QC and get alignment stats (QC)

3. Quantify

  * **FeatureCounts**: count reads per feature

4. Report

  * calculate TPM values
  * generate report containing count data (for use in downstream DGE analysis e.g. limma/voom)


<a name="todo"/>

## To Do

* unit tests and example data

* estimate library complexity for each sample (try using [preseq](http://smithlabresearch.org/software/preseq/))

* consolidate calls to MultiQC to generate only one report (reports are currently generated after both QC and alignment)

* incorporate existing standalone scripts to:

	* autogeneration of (BigsWigs and) [JBrowse](http://jbrowse.org/) instances

	* setup local [Degust](http://degust.erc.monash.edu/)

	* generate heatmaps to visualise sample similarity (for identifying batch effects)

* fully integrate DupRadar and coverage checks (these aren't reported in the REPORT files and are checked manually)


<a name="installation"/>

## Installation

* Use Anaconda to manage pipeline dependencies (*optional*)

```
conda create -n python2 python=2.7 anaconda

source activate python2
```

* Download the repo

```
git clone https://github.com/will-rowe/rnaseq && cd rnaseq
```

* Edit GLOBALS.py file to edit setup variables

```
vim ./scripts/GLOBALS.py
```

* Check (and install) required software

```
more requirements.txt
```

* Edit the dupRadar wrapper script (*optional*)

```
vim ./scripts/dupRadar.sh
```

* Run the setup.py (use the --user option if you are not in a virtual environment or have no root privileges)

```
python setup.py install --user
```

* Try calling the program

```
rnaseq --help
```


<a name="running"/>

## Running

The basic command to run the pipeline is:

```
rnaseq [OPTIONS] reference file_list
```

**reference**:	reference to be used for read alignment / annotation (must be in GFF3 format and include the fasta sequence)

**file_list**:	the list of fastq files to be processed (must be newline separated text file containing full paths to files)


<a name="output"/>

## Output

A report is generated for each input file. To combine the count data from each report into a single *tsv* file (to use in R for DGE / upload to Degust), try something like:

1. remove comment lines from each file and get the feature names

```
for i in *.REPORT; do sed -i -e '1,5d' $i; done

cut -f1 filename.REPORT > features.tmp && sed -i "1s/^/FEATURE\n/" features.tmp
```

2. join the 4th columns for a combined count file - making sure to add in the file name to the top of the column

```
for i in *.REPORT; do cut -f4 $i > $i.countsonly; sed -i "1s/^/$i\n/" $i.countsonly; done

paste features.tmp *.countsonly -d '\t' > D23-counts.tsv
```

3. join the 5th columns for a combined TPM file

```
for i in *.REPORT; do cut -f5 $i > $i.tpmsonly; sed -i "1s/^/$i\n/" $i.tpmsonly; done

paste features.tmp *.tpmsonly -d '\t' > D23-tpms.tsv
```


<a name="issues"/>

## Known Issues

* The pipeline checks for required software but does not check all version numbers - this needs to be checked by the user (using the requirements.txt file)

* A failed trimmomatic run won't be reported in the log

* No check for Kraken database included at pipeline initialisation

* Issues when installing Kraken on MacOS - try referring to [this](https://groups.google.com/forum/#!msg/kraken-users/MMBbOSPVeQM/8EyqRrvaAx8J;context-place=msg/kraken-users/0i2m81J_1uw/_zk1zUB3npsJ)

* Issue when installing rnaseq on Ubuntu 16:

```
freetype: no  	[The C/C++ header for freetype2 (ft2build.h)
            	could not be found.  You may need to install the
            	development package.]
```

you might need to install pkg-config:
```
conda install --name python2 pkg-config
```
