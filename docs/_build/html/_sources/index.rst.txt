.. MetaboDirect documentation master file, created by
   sphinx-quickstart on Thu Jun 17 14:26:12 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================================
Welcome to MetaboDirect's documentation!
========================================

|pypi_badge|

.. |pypi_badge| image:: https://img.shields.io/pypi/v/metabodirect?style=plastic
    :target: https://pypi.org/project/metabodirect/
    
------------
Introduction
------------

**MetaboDirect** is a Python and R based pipeline for the analysis of Direct Injection FT-ICR Mass Spectrometry data.
The MetaboDirect pipeline takes a *report file*, generated using Formularity (Tolić et al., 2017), and a *sample information file* as inputs and automatically performs all the analysis described below including: sample filtering, m/z filtering, normalization of intensities, thermodynamic index calculation, annotation of molecular formulas using the KEGG database (Kanehisa & Goto, 2000), statistical analysis and construction of transformation networks using Cytoscape (Shannon et al., 2003).


.. toctree::
   :maxdepth: 2
   
   installation
   quickstart
   pipeline
   guide


------------------
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
