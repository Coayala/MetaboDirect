# User Guide for MetaboDirect

## Introduction

**MetaboDirect** is a Python and R based pipeline for the analysis of Direct Injection FT-ICR Mass Spectrometry data.
The MetaboDirect pipeline takes a Formularity (Tolić et al., 2017) *report* and a *sample information* file as inputs and automatically performs all the analysis described below including: sample filtering, m/z filtering, normalization of intensities, thermodynamic index calculation, annotation of molecular formulas using the KEGG database (Kanehisa & Goto, 2000), statistical analysis and construction of transformation networks using Cytoscape (Shannon et al., 2003)

## Quick Start

To quickly run the MetaboDirect pipeline use the following command. Read below for a more detailed information about the *input data* and the analysis options that are offered.

```
/path/to/MetaboDirect/run_metabodirect.py DATA_FILE METADATA_FILE -g GROUPING_VARIABLE
```
Information about the arguments can be obtained using the -h/--help function.

```
/path/to/MetaboDirect/run_metabodirect.py -h
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

The **MetaboDirect** pipeline includes 5 major steps: data pre-processing, data diagnostics, data exploration, statistical analysis, and transformation network analysis. In addition an optional script `normalization_test.py` can be run before **MetaboDirect** to help the user choose the best normalization method for the data.

#### Input data

The input data for **MetaboDirect** is the a .csv report generated by default by Formularity. If Formularity was not used, the data file should be arranged to have columns with the **exact same names** that are shown below. 
*Mass1, Mass2, ...* refer to the m/z values detected by the software (i.e. each peak) while *Sample1, Sample2, ..* refer to the intesity of each peak for each of the samples.
An example dataset is included in the `/path/to/MetaboDirect/data/` directory with the name `Report.csv`.

|Mass|C|H|O|N|C13|S|P|Na|El_comp|Class|NeutralMass|Error_ppm|Candidates|*Sample1*|*Sample2*|*Sample3*|*...*|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|*Mass1*|4|4|4|4|0|0|0|0||NA|111.634719|0|NA|6.237|0|0|
|*Mass2*|5|2|2|0|0|0|0|0||NA|111.712035|0|NA|0|6.343|6.166|
|*Mass3*|3|6|2|2|0|1|0|0||NA|112.125136|0|NA|7.549|7.363|6.75|
|*Mass4*|5|6|3|0|0|0|1|0||NA|112.3957945|0|NA|0|0|6.145|
|*Mass5*|6|2|3|0|0|0|1|0||NA|112.457043|0|NA|0|6.133|0|
|*...*|||||||||||||||||

#### Sample information file

The sample informatio file is a .csv file that has one column called *SampleID* with the names of all of the samples that are present in the report file. Please make sure that the sample names in the *input data* and the *sample information file* are **exactly the same**. At least one other column must be present in the sample information file and must contain information used to group the data for plotting and for the statistical analysis. Multiple grouping variables can be present in this file but only two can be used simultaneously in **MetaboDirect**. When running the pipelines the grouping variables can be defined with the `-g` option using the **exact name** that it is on this file.  Additionally, please use only letters (Aa-Zz), numbers (0-9) and underscores (\_) for both the **sample names** and the **grouping variables**.
An example file is included in the `/path/to/MetaboDirect/data/` directory with the name `metadata.csv`.

|SampleID|Grouping_var1|Grouping_var2|Grouping_var3
|---|---|---|---|
|*Sample1*|A|M|X|
|*Sample2*|A|N|Y|
|*Sample3*|B|M|Y|
|*Sample4*|B|N|X|

### (Optional Step) Test normalization methods

A companion script called `normalization_test.py` is included to help in the decision of which normalization method to use. This script uses a procedure the Statistical Procedure for the Analysis of Normalization Strategies (SPANS) (Webb‐Robertson et al., 2011), which has been previously demonstrated to work well with FT-ICR MS data (Thompson et al., 2021).
Information about the arguments can be obtained using the -h/--help function.

```
/path/to/MetaboDirect/normalization_test.py -h
```

### 1. Data pre-processing

![Data pre-processing](/docs/pre_processing.png)

### 2. Data diagnostics

![Data diagnostics](/docs/diagnostics.png)

### 3. Data exploration

![Data exploration](/docs/exploration.png)

### 4. Statistical Analysis

![Statistical Analysis](/docs/statistics.png)

### 5. Transformation Networks

![Transformation Networks](/docs/transformations.png)


## References

- Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Research, 28(1), 27-30. https://doi.org/10.1093/nar/28.1.27 
- Shannon, P., Markiel, A., Ozier, O., Baliga, N. S., Wang, J. T., Ramage, D., Amin, N., Schwikowski, B., & Ideker, T. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome research, 13(11), 2498-2504. https://doi.org/10.1101/gr.1239303 
- Thompson, A. M., Stratton, K. G., Bramer, L. M., Zavoshy, N. S., & McCue, L. A. (2021). Fourier transform ion cyclotron resonance mass spectrometry (FT-ICR-MS) peak intensity normalization for complex mixture analyses [https://doi.org/10.1002/rcm.9068]. Rapid Communications in Mass Spectrometry, 35(9), e9068. https://doi.org/https://doi.org/10.1002/rcm.9068
- Tolić, N., Liu, Y., Liyu, A., Shen, Y., Tfaily, M. M., Kujawinski, E. B., Longnecker, K., Kuo, L.-J., Robinson, E. W., Paša-Tolić, L., & Hess, N. J. (2017). Formularity: Software for Automated Formula Assignment of Natural and Other Organic Matter from Ultrahigh-Resolution Mass Spectra. Analytical Chemistry, 89(23), 12659-12665. https://doi.org/10.1021/acs.analchem.7b03318
- Webb‐Robertson, B. J. M., Matzke, M. M., Jacobs, J. M., Pounds, J. G., & Waters, K. M. (2011). A statistical selection strategy for normalization procedures in LC‐MS proteomics experiments through dataset‐dependent ranking of normalization scaling factors. Proteomics, 11(24), 4736-4741.


