# RNA-seq pipeline
This is the repo for the Hinton Lab RNA-seq pipeline

-----------


## General updates

This pipeline is still under active development --- extensive readme and tests to follow


## Contact
### Will Rowe

* Website: [Will Rowe](https://github.com/will-rowe)
* Lab website: [hintonlab.com](http://www.hintonlab.com)
* Email: [Will Rowe @ Liverpool](will.rowe@liverpool.ac.uk)


## Installation

1. Download the repo

`git clone https://github.com/will-rowe/rnaseq && cd rnaseq`


2. Edit GLOBALS.py file to edit setup variables

`vim ./scripts/GLOBALS.py`


3. Edit dupRadar script if needed

`vim ./scripts/dupRadar.sh`


4. Run the setup.py


`python setup.py install --user`


5. Call the program

`rnaseq --help`
