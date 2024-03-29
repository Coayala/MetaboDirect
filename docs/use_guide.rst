============
User's Guide
============

|pypi_badge|

.. |pypi_badge| image:: https://img.shields.io/pypi/v/metabodirect?style=plastic
    :target: https://pypi.org/project/metabodirect/

.. _metabod:

----------------------
**metabodirect**
----------------------

**MetaboDirect** requires two files, a report from Formularity (data file) and a sample information file (metadata file), and at least one grouping variable (defined in the metadata file).

Information about the arguments can be obtained using the ``-h/--help`` function. A detailed description of each argument is provided below.

.. code-block:: none

	metabodirect -h

.. code-block:: none

	usage: metabodirect [-h] -g STR [STR ...] [-o OUTDIR] [-m FLOAT FLOAT] [-p INT] [-e FLOAT] [-f STR STR] [-b STR] [-k]
			    [-v] [-n STR] [--norm_subset STR] [--subset_parameter FLOAT] [--log_transform] [-t] [-c]
			    DATA METADATA

	Program for running all the MetaboDirect analysis pipeline

	positional arguments:
	  DATA                  Name of the file with the Direct Injection MS data in .csv format
	  METADATA              Name of the file with the sample information (metadata) in .csv format

	optional arguments:
	  -h, --help            show this help message and exit
	  -g STR [STR ...], --group STR [STR ...]
				Grouping variables for coloring and faceting figures (Max 2) (default: None)
	  -o OUTDIR, --outdir OUTDIR
				Output directory (default: MetaboDirect_output)
	  -m FLOAT FLOAT, --mass_filter FLOAT FLOAT
				Range to filter m/z data (min_mz, max_mz). The pipeline will not filter m/z values by default
				(default: None)
	  -p INT, --peak_filter INT
				Minimum number of samples a peak must be present to be conserved for the analysis (default: 2)
	  -e FLOAT, --error_filter FLOAT
				Max error (e) allowed in formula assignment. Peaks with |error| > e will be removed from the
				analysis (default: 0.5)
	  -f STR STR, --filter_by STR STR
				Filter samples based on metadata. First enter the name of the feature, followed by the values
				associated with the samples you want to keep in the analysis.(Example -f Habitat Bog,Palsa)
				(default: None)
	  -b STR, --biochem_key STR
				File with the biochemical key to use for the transformation network (default: Default key)
	  -k, --kegg_annotation
				Set this option to perform annotation of the molecular formulas usingthe KEGG database
				(default: False)
	  -v, --version         show program's version number and exit

	Normalization methods:
	  Options to define how data normalization will be carried out

	  -n STR, --norm_method STR
				Available methods to normalize data are: 'mean', 'median', 'zscore', 'sum', 'max', 'minmax',
				'binary', 'none' (default: max)
	  --norm_subset STR     Subset of the data to use for normalization purpouses. Available subset methods: ALL, LOS,
				PPP. LOS uses peaks in the top L order statistics, PPP uses peaks having a minimum percentage
				of observed values. (default: ALL)
	  --subset_parameter FLOAT
				If using a sample subset for nomalization, this parameter defines the subsample of peaks that
				will be used for normalization. If not defined, the default values will be 0.3 for LOS and 0.5
				for PPP (default: None)
	  --log_transform       Set this option to log transform the data. (Program will fail if there are peaks with
				intensities of 0. Consider tranforming this values into 1 if log transformation is desired
				(default: False)

	Transformation network options:
	  Options to control wheter transformations will be calculated and if networks will be constructed

	  -t, --calculate_transformations
				Set this option to calculate transformations based on biochemical key (default: False)
	  -c, --create_networks
				Set this option to build transformation networks based on transfomations calculatedwith the
				biochemical key (this options turns -t automatically) (default: False)

++++++++++++++++++++++++++
Input data file (``DATA``)
++++++++++++++++++++++++++

The input data for **MetaboDirect** is the .csv report file generated by default by Formularity. If Formularity was not used, any data file can be arranged to have columns with the **exact same names** that are shown in the table below. 

.. csv-table::
	:header: "Mass", "C", "H", "O", "N", "C13", "S", "P", "Na", "El_comp", "Class", "NeutralMass", "Error_ppm", "Candidates", "*Sample1*", "*Sample2*", "*Sample3*", "*...*"

	"*Mass1*", "4", "4", "4", "4", "0", "0", "0", "0", "", "NA", "111.634719", "0", "NA", "6.237", "0", "0"
	"*Mass2*", "5", "2", "2", "0", "0", "0", "0", "0", "", "NA", "111.712035", "0", "NA", "0", "6.343", "6.166"
	"*Mass3*", "3", "6", "2", "2", "0", "1", "0", "0", "", "NA", "112.125136", "0", "NA", "7.549", "7.363", "6.75"
	"*Mass4*", "5", "6", "3", "0", "0", "0", "1", "0", "", "NA", "112.3957945", "0", "NA", "0", "0", "6.145"
	"*Mass5*", "6", "2", "3", "0", "0", "0", "1", "0", "", "NA", "112.457043", "0", "NA", "0", "6.133", "0"
	"*...*", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""

*Mass1, Mass2, ...* refer to the m/z values detected by the software (i.e. each peak) while *Sample1, Sample2, ..* refer to the name of each sample.
An example dataset is included in the MetaboDirect repository ``example`` `directory <https://github.com/Coayala/MetaboDirect/tree/main/example>`_.

++++++++++++++++++++++++++++++++++++++
Sample information file (``METADATA``)
++++++++++++++++++++++++++++++++++++++

The sample information file (or metadata file) is a .csv file that has one column called *SampleID* with the names of all of the samples that are present in the report file. Please make sure that the sample names in the *input data* and the *sample information file* are **exactly the same**. At least one other column must be present in the sample information file and must contain information used to group the data for plotting and for the statistical analysis. Multiple grouping variables can be present in this file but only two can be used simultaneously in **MetaboDirect**. When running the pipelines the grouping variables can be defined with the `-g` option using the **exact name** that it is on this file.  Additionally, please use only letters (Aa-Zz), numbers (0-9) and underscores ( \_ ) for both the **sample names** and the **grouping variables**.
An example file is included in the ``example`` `folder <https://github.com/Coayala/MetaboDirect/tree/main/example>`_ with the name `metadata.csv`.

.. csv-table::
	:header: "SampleID", "Grouping_var1", "Grouping_var2", "Grouping_var3"
	
	"*Sample1*", "A", "M", "X"
	"*Sample2*", "A", "N", "Y"
	"*Sample3*", "B", "M", "Y"
	"*Sample4*", "B", "N", "X"
	"*Sample5*", "A", "N", "Z"
	
++++++++++++++++++++++++++++++++++++++++
Output directory (``-o`` | ``--outdir``)
++++++++++++++++++++++++++++++++++++++++

The name of directory where all the generated plots, tables and scripts will be saved. If it is not defined the directory will be named MetaboDirect_output by default.

++++++++++++++++++++++++++++++++++++++++
Grouping variable (``-g`` | ``--group``)
++++++++++++++++++++++++++++++++++++++++

This option accepts up to two grouping variables (e.g. ``-g Grouping_var1`` or ``-g Grouping_var1 Grouping_var2``) whose names are **exactly the same** as they appear in the columns of the metadafile. The first grouping variable will be used for giving colors to the plots generated. Both variables will be used for the statistical analysis and the pairwise comparisons.

+++++++++++++++++++++++++++++++++++++++++
Filter samples (``-f`` | ``--filter_by``)
+++++++++++++++++++++++++++++++++++++++++

This option takes two arguments: **1)** a variable from the metadata file and **2)** values from that variable column that we want to keep in the analysis. For example ``-f Grouping_var3 X``, will keep just the samples for whom the Groupin_var3 is equal to "X". Multiple values for the same variable can be defined separated by commas (without spaces) (i.e. ``-g Grouping_var3 X,Z``).

++++++++++++++++++++++++++++++++++++++++
Mass filter (``-m`` | ``--mass_filter``)
++++++++++++++++++++++++++++++++++++++++

This option takes two arguments: lower and an upper m/z limits. Peaks with m/z (masses) outside of its limits will be filtered out and not considered in the analysis.

++++++++++++++++++++++++++++++++++++++++
Peak filter (``-e`` | ``--error_filter``)
++++++++++++++++++++++++++++++++++++++++

This option is to determine the maximum error that is allowed from formula assignment.

++++++++++++++++++++++++++++++++++++++++
Error filter (``-p`` | ``--peak_filter``)
++++++++++++++++++++++++++++++++++++++++

This option is for specified the minimum number of samples a peak must be present to be conserved for the analysis.

+++++++++++++++++++++++++++++++++++++++++++++++++
Normalization method (``-n`` | ``--norm_method``)
+++++++++++++++++++++++++++++++++++++++++++++++++

This option defines which normalization method will be used to normalize the intensities (*I*). It can take one of the following options for *i* samples and *j* peaks.
Normalization methods are based on the ones used by Kitson, et al. (2021) and Thompson, et al. (2021):

.. csv-table::
	:header: "Normalization method", "Formula"
	
	"``max``", ":math:`NormIntensity_{i,j} = \frac{I_{i,j}}{max(I)_{i}}`"
	"``minmax``", ":math:`NormIntensity_{i,j} = \frac{I_{i,j} - min(I)_{i}}{max(I)_{i} - min(I)_{i}}`"
	"``mean``", ":math:`NormIntensity_{i,j} = \frac{I_{i,j} - mean(I)_{i}}{max(I)_{i} - min(I)_{i}}`"
	"``median``", ":math:`NormIntensity_{i,j} = \frac{I_{i,j} - median(I)_{i}}{max(I)_{i} - min(I)_{i}}`"
	"``sum``", ":math:`NormIntensity_{i,j} = \frac{I_{i,j}}{\sum{I}_{i}}`"
	"``zscore``", ":math:`NormIntensity_{i,j} = \frac{I_{i,j} - mean(I)_{i}}{std.dev(I)_{i}}`"
	"``none``", ":math:`NormIntensity_{i,j} = InputData_{i,j}`"

+++++++++++++++++++++++++++++++++++++++++++++++
Normalization subset method (``--norm_subset``)
+++++++++++++++++++++++++++++++++++++++++++++++

If a normalization method other than ``binary`` or ``none`` is selected it is possible to use only a fraction of the peaks to calculate the normalization factors (normalization will still be applied to all the dataset). Possible subset methods are :

.. csv-table::
	:header: "Subset method", "Description"
	
	"``ALL``", "Use all present peaks to calculate normalization factors"
	"``LOS``", "Use a percentage of peaks in the top L order statistics"
	"``PPP``", "Uses peaks that are present in more than minimum percentage of samples"
	
The option ``--subset_parameter`` defines the percentage of peaks that will be used in ``LOS`` or the minimum percentage of samples that a peak must be present for ``PPP``.

++++++++++++++++++++++++++++++++++++
Subset parameter (``--norm_subset``)
++++++++++++++++++++++++++++++++++++

This option is only needed when ``LOS`` of ``PPP`` are selected as normalization methods. It defines either the minimum percentage of samples a peaks need to bre present to be considered (``PPP``) or the percentage of top peaks that will be used (``LOS``).

++++++++++++++++++++++++++++++++++++++++++++++++
KEGG annotation (``-k`` | ``--kegg_annotation``)
++++++++++++++++++++++++++++++++++++++++++++++++

This is an optional step as it may take a long time (~ couple of hours) depending on the number of peaks present in the data. If this option is present, peaks will be annotated with the KEGG database (Pathway, Module, Brite, etc.) based on their molecular formula.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Calculate transformations (``-t`` | ``--calculate_transformations``)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This option define whether or not a molecular transformations between the peaks will be calculated based on their mass differences. If this option is selected, **MetaboDirect** will end after generating the transformation files. Transformation files will be located in ``./$outdir/6_transformations/transf_by_sample``.

++++++++++++++++++++++++++++++++++++++++++++++++++++++
Create networks (``-c`` | ``--create_networks``)
++++++++++++++++++++++++++++++++++++++++++++++++++++++

If this option is selected, it will automatically turn on the option ``-t``. After the transformation files are generated, transformation networks will be built. This step requires Cytoscape (version 3.8 and above) to be installed in the machine. **MetaboDirect** will ask the user to open Cytoscape when required in order to construct the networks. When prompted in the screen, please open Cytoscape and then hit enter to continue with the analysis.
	

.. _testnorm:

----------------------
**test_normalization**
----------------------

This a companion script that can be used to help choosing the best normalization method for the data using the SPANS method.

Information about the arguments can be obtained using the ``-h/--help`` function. A detailed description of each argument is provided below.

.. code-block:: none

	test_normalization -h

.. code-block:: none

	usage: test_normalization [-h] [-f STR STR] [--log_transform] DATA METADATA GROUP

	Program for running all the MetaboDirect analysis pipeline

	positional arguments:
	DATA                  Name of the file with the DI-MS data in .csv format
	METADATA              Name of the file with the sample information (metadata) in tabular format
	GROUP                 Grouping variables to test for normalization significance

	optional arguments:
	-h, --help            show this help message and exit
	-f STR STR, --filter_by STR STR
							Filter samples based on metadata. First enter the name of the feature,followed by the values
							associated with the samples you want to keep in the analysis.(Example -f Habitat Bog,Palsa)
							(default: None)
	--log_transform       Set this if you plan to log transform your data before normalization (default: False)

++++++++++++++++++++++++++
Input data file (``DATA``)
++++++++++++++++++++++++++

The same input data that will be used for **MetaboDirect**. A .csv report file generated by default by Formularity. If Formularity was not used, any data file can be arranged to have columns with the **exact same names** that are shown above for :ref:`metabod`. 

++++++++++++++++++++++++++++++++++++++
Sample information file (``METADATA``)
++++++++++++++++++++++++++++++++++++++

The same input data that will be used for **MetaboDirect**. A .csv file that has one column called *SampleID* with the names of all of the samples that are present in the report file. Please make sure that the sample names in the *input data* and the *sample information file* are **exactly the same**.

++++++++++++++++++++++++++++++++++++++++
Grouping variable (``GROUP``)
++++++++++++++++++++++++++++++++++++++++

The grouping variable that will be tested for significance . It names should be **exactly the same** as they appear in the columns of the metadafile.

+++++++++++++++++++++++++++++++++++++++++
Filter samples (``-f`` | ``--filter_by``)
+++++++++++++++++++++++++++++++++++++++++

This option takes two arguments: **1)** a variable from the metadata file and **2)** values from that variable column that we want to keep in the analysis. For example ``-f Grouping_var3 X``, will keep just the samples for whom the Groupin_var3 is equal to "X". Multiple values for the same variable can be defined separated by commas (without spaces) (i.e. ``-g Grouping_var3 X,Z``).

++++++++++++++++++++++++++++++++++++++++
Log transform data (``--log_transform``)
++++++++++++++++++++++++++++++++++++++++

If this option is used, data will be log transformed before testing for significance.

.. _createnet:

-------------------
**create_networks**
-------------------

This a companion script that can be used to build the transformation networks using the transformation files created with the option ``-t``.

Information about the arguments can be obtained using the ``-h/--help`` function. A detailed description of each argument is provided below.

.. code-block:: none

	create_networks -h

.. code-block:: none

	usage: create_networks [-h] OUTDIR METADATA STR [STR ...]

	Program for creating molecular transformation networks, based on previously calculated transformations

	positional arguments:
	OUTDIR      Output directory used to create networks with metabodirect and the -t option
	METADATA    Metadata file used in the analysis, if a filtered metadata was generated please enter that one
	GROUP         Grouping variables for coloring and faceting figures (Max 2)

	optional arguments:
	-h, --help  show this help message and exit

+++++++++++++++++++++++++++++
Output directory (``OUTDIR``)
+++++++++++++++++++++++++++++

It needs to be the same output directory that was used during the original run of the **MetaboDirect** pipeline with the ``-t`` option.

++++++++++++++++++++++++++++++++++++++
Sample information file (``METADATA``)
++++++++++++++++++++++++++++++++++++++

The same metadata file that was used during the original run of the **MetaboDirect** pipeline with the ``-t`` option.

++++++++++++++++++++++++++++++++++++++++
Grouping variable (``GROUP``)
++++++++++++++++++++++++++++++++++++++++

The grouping variable that will be used to compare the network statistics. It names should be **exactly the same** as they appear in the columns of the metadafile.


++++++++++
References
++++++++++

- Kitson, E., Kew, W., Ding, W., & Bell, N. G. A. (2021). PyKrev: A Python Library for the Analysis of Complex Mixture FT-MS Data. Journal of the American Society for Mass Spectrometry, 32(5), 1263-1267. https://doi.org/10.1021/jasms.1c00064 
- Thompson, A. M., Stratton, K. G., Bramer, L. M., Zavoshy, N. S., & McCue, L. A. (2021). Fourier transform ion cyclotron resonance mass spectrometry (FT-ICR-MS) peak intensity normalization for complex mixture analyses [https://doi.org/10.1002/rcm.9068]. Rapid Communications in Mass Spectrometry, 35(9), e9068. https://doi.org/https://doi.org/10.1002/rcm.9068
