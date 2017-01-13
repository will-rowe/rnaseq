# RNA-seq pipeline
This is the repo for the Hinton Lab RNA-seq pipeline

-----------


## General updates

This repo only contains the pipeline wrapper and the QC / GFF modules - more to follow once testing is complete


## Contact
### Will Rowe

* Website: [Will Rowe](https://github.com/will-rowe)
* Lab website: [hintonlab.com](http://www.hintonlab.com)
* Email: [Will Rowe @ Liverpool](will.rowe@liverpool.ac.uk)
* Twitter: [@will__rowe](https://twitter.com/will__rowe)


## Installation

1. Download the repo

`git clone https://github.com/will-rowe/rnaseq && cd rnaseq`

2. Edit `./scripts/GLOBALS.py` to edit setup variables

3. Edit `./scripts/dupRadar.sh` if needed

4. Run the setup.py

`python setup.py install --user`

5. Call the program

`rnaseq --help`
