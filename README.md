# StPeter
MS2 spectral intensity-based label-free quantification of LC-MS/MS data

NOTE: This is an early prototype of the software. The current version is compiled into the [Trans-Proteomic Pipeline](https://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20%28TPP%29/)

## Installation
`> virtualenv env`

`> source env/bin/activate`

`> pip install -r requirements.txt`

## Usage
Standard 1% FDR cutoff and 0.4 Da MS2 tolerance

`> python StPeter.py myfile.prot.xml`

To see additional options:

`> python StPeter.py -h`
