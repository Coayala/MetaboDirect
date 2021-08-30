============
Installation
============

|pypi_badge|

.. |pypi_badge| image:: https://img.shields.io/pypi/v/metabodirect?style=plastic
    :target: https://pypi.org/project/metabodirect/

MetaboDirect can be installed directly from `PyPi <https://pypi.org/project/metabodirect/0.1.1/>`_ using:

.. code-block:: none

	pip install metabodirect


Additionally it can be installed from source by cloning its `GitHub repository <https://github.com/Coayala/MetaboDirect>`_:

.. code-block:: bash
	
	git clone https://github.com/Coayala/MetaboDirect.git
	cd MetaboDirect
	python setup.py install


Required modules
----------------

MetaboDirect requires Python (3.5 and above), R (4 and above) and Cytoscape (3.8 and above) with the following libraries/modules:

Python
++++++

- argparse
- numpy
- pandas
- seaborn
- more-itertools
- py4cytoscape

R
++

- tidyverse
- RColorBrewer
- vegan
- ggnewscale
- ggpubr
- vegan
- KEGGREST
- factoextra
- UpSetR
- pmartR (for normalization tests)
- SYNCSA

Cytoscape
+++++++++

- FileTransfer