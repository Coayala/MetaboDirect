# MetaboDirect for agentic AI

This branch is intended for a version of MetaboDirect that can be used with an AI agent

It needs to be installed from source by cloning the `agent` branch from the [GitHub repository](https://github.com/Coayala/MetaboDirect/tree/agent)

```
git clone -b agent https://github.com/Coayala/MetaboDirect.git
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
- KEGGREST
- factoextra
- UpSetR
- pmartR (for normalization tests)
- SYNCSA
- ggvenn
- ggrepel
- glue
- jsonlite

### Cytoscape

- FileTransfer

## Usage

Information about the arguments can be obtaining using the option -h/--help

```
metabodirect -h
```
For more information please check the [User Manual](https://metabodirect.readthedocs.io/en/latest/index.html#).


