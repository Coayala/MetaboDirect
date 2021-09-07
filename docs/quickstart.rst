===========
Quick start
===========

|pypi_badge|

.. |pypi_badge| image:: https://img.shields.io/pypi/v/metabodirect?style=plastic
    :target: https://pypi.org/project/metabodirect/

To quickly run the MetaboDirect pipeline use the following command. 

.. code-block:: none

	metabodirect DATA_FILE METADATA_FILE -g GROUPING_VARIABLE

Check :ref:`metabod` in the User's Guide to learn more about the required *input data* and the analysis options that are offered.

Information about the arguments can be obtained using the -h/--help function.

.. code-block:: none

	metabodirect -h

-----------------------
Example using test data
-----------------------

Example data can be downloaded from the MetaboDirect repository ``example`` `directory <https://github.com/Coayala/MetaboDirect/tree/main/example>`_, or from the command line with:

.. code-block:: none

	# Report file
	wget https://raw.githubusercontent.com/Coayala/MetaboDirect/main/example/Report.csv --no-check-certificate

	# Metadata file
	wget https://raw.githubusercontent.com/Coayala/MetaboDirect/main/example/metadata.csv --no-check-certificate


Try analizing example data using:

.. code-block:: none

	metabodirect Report.csv metadata.csv -o test -m 200 400 -g Habitat Depth -t -k  
