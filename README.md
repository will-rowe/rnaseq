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
