# MetaboDirect

[![DOI](https://zenodo.org/badge/356707576.svg)](https://zenodo.org/badge/latestdoi/356707576)

A comprehensive command-line based pipeline for the analysis of direct injection FT-ICR mass spectrometry data.

## Installation

MetaboDirect can be installed directly from [PyPi](https://pypi.org/project/metabodirect/0.1.1/) using:

```
pip install metabodirect
```

Additionally it can be installed from source by cloning its [GitHub repository](https://github.com/Coayala/MetaboDirect)

```
git clone https://github.com/Coayala/MetaboDirect.git
cd MetaboDirect
python setup.py install
```

MetaboDirect requires Python (3.5 and above), R (4 and above) and Cytoscape (3.8 and above) with the following libraries/modules:

### Python

- argparse
- numpy
- pandas
- seaborn
- more-itertools
- py4cytoscape
- statsmodels

### R

- tidyverse
- RColorBrewer
- vegan
- ggnewscale
- ggpubr
- ggvenn
- KEGGREST
- factoextra
- UpSetR
- pmartR (for normalization tests)
- SYNCSA
- ggvenn

### Cytoscape

- FileTransfer

## Usage

Information about the arguments can be obtaining using the option -h/--help

```
metabodirect -h
```
For more information please check the [User Manual](https://metabodirect.readthedocs.io/en/latest/index.html#).


