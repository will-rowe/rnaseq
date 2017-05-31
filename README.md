# RNA-seq pipeline
This is the repo for the Hinton Lab RNA-seq pipeline

-----------


## General updates

This pipeline is still under active development --- extensive readme and tests to follow


## Contact
### Will Rowe

* Website: [Will Rowe](https://will-rowe.github.io)
* Lab website: [hintonlab.com](http://www.hintonlab.com)
* Email: [Will Rowe @ Liverpool](will.rowe@liverpool.ac.uk)


## Installation

Optional - We use Anaconda to manage pipeline dependencies:

`conda create -n python2 python=2.7 anaconda`

`source activate python2`

1. Download the repo

`git clone https://github.com/will-rowe/rnaseq && cd rnaseq`


2. Edit GLOBALS.py file to edit setup variables

`vim ./scripts/GLOBALS.py`


3. Edit dupRadar script if needed

`vim ./scripts/dupRadar.sh`


4. Run the setup.py (use the --user option if you are not in a virtual environment or have no root privileges)

`python setup.py install --user`


5. Call the program

`rnaseq --help`



## Output

A report is generated for each input file. To combine these into a single tsv, try something like:

1. remove comment lines from each file and get the feature names
for i in *.REPORT; do sed -i -e '1,5d' $i; done
cut -f1 filename.REPORT > features.tmp && sed -i "1s/^/FEATURE\n/" features.tmp

2. join the 4th columns for a combined count file - making sure to add in the file name to the top of the column
for i in *.REPORT; do cut -f4 $i > $i.countsonly; sed -i "1s/^/$i\n/" $i.countsonly; done
paste features.tmp *.countsonly -d '\t' > D23-counts.tsv


3. join the 5th columns for a combined TPM file
for i in *.REPORT; do cut -f5 $i > $i.tpmsonly; sed -i "1s/^/$i\n/" $i.tpmsonly; done
paste features.tmp *.tpmsonly -d '\t' > D23-tpms.tsv
