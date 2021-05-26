# User Guide for MetaboDirect
## Introduction

MetaboDirect is a Python and R based pipeline for the analysis of Direct Injection FT-ICR Mass Spectrometry data.
The MetaboDirect pipeline takes a Formularity (Tolić et al., 2017) *report* and a *sample information* file as inputs and automatically performs all the analysis described below including: sample filtering, m/z filtering, normalization of intensities, thermodynamic index calculation, annotation of molecular formulas using the KEGG database (Kanehisa & Goto, 2000), statistical analysis and construction of transformation networks using Cytoscape (Shannon et al., 2003)

## Quick Start

To run the MetaboDirect pipeline use the following command

```
$METABODIRECT_HOME/run_metabodirect.py DATA_FILE METADATA_FILE -g GROUPING_VARIABLE
```
Information about the arguments can be optained using the -h/--help function

```
$METABODIRECT_HOME/run_metabodirect.py -h
```

```
usage: run_metabodirect.py [-h] [-o OUTDIR] [-m INT INT] -g STR [STR ...] [-f STR STR] [-b STR] [-t] [-k] [-n STR] [--norm_subset STR] [--subset_parameter STR] [--log_transform]
                           DATA METADATA

Program for running all the MetaboDirect analysis pipeline

positional arguments:
  DATA                  Name of the file with the DI-MS data in .csv format
  METADATA              Name of the file with the sample information (metadata) in .csv format

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Output directory (default: MetaboDirect_output)
  -m INT INT, --mass_filter INT INT
                        Range to filter m/z data (min_mz, max_mz). The pipeline will not filter m/z values by default (default: None)
  -g STR [STR ...], --group STR [STR ...]
                        Grouping variables for coloring and faceting figures (Max 2) (default: None)
  -f STR STR, --filter_by STR STR
                        Filter samples based on metadata. First enter the name of the feature,followed by the values associated with the samples you want to keep in the
                        analysis.(Example -f Habitat Bog,Palsa) (default: None)
  -b STR, --biochem_key STR
                        File with the biochemical key to use for the transformation network (default: Default key)
  -t, --transformation_analysis
                        Set this option to perform a transformation netwokr analysis of the samples (default: False)
  -k, --kegg_annotation
                        Set this option to perform annotation of the molecular formulas usingthe KEGG database (default: False)

Normalization methods:
  Options to define how data normalization will be carried out

  -n STR, --norm_method STR
                        Available methods to normalize data (default: max)
  --norm_subset STR     Subset of the data to use for normalization purpouses. LOS uses peaks in the top L order statistics, PPP uses peaks having a minimum percentage of observed
                        values. (default: ALL)
  --subset_parameter STR
                        If using a sample subset for nomalization, this parameter defines the subsample of peaks that will be used for normalization.If not defined, the default values
                        will be 0.3 for LOS and 0.5 for PPP (default: None)
  --log_transform       Set this option to log transform the data. (Program will fail if there are peaks with intensities of 0. Consider tranforming this values into 1 if log
                        transformation is desired (default: False)
```

## MetaboDirect pipeline

The MetaboDirect pipeline includes 5 major steps: data pre-processing, data diagnostics, data exploration, data statistics and transformation network analysis.



## (Optional)


## References


