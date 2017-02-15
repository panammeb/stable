DATASETS
-----------------------------------------------------------------------
Testing PANAM

In order to test PANAM, a full 454Roche format dataset is provided.
This dataset contains 4 files:
- a sequence file : 454roche.seq
- a quality file : 454roche.qual
- a barcode file : 454roche.bar
- a configuration file : panam.ini

For these who would test PANAM with MiSeq format data, only a configuration
file is provided, as an example. Due to space quota limitations on github
hosting, it is not possible to provide MiSeq format input data files.

-----------------------------------------------------------------------
How to set datasets to test PANAM ?

Configuration files are built to allow you testing PANAM assuming that 
you have installed PANAM in a subdirectory named panam in your home
directory, and that datasets have been deployed in a subdirectory named
454roche (or miseq, if applicable) in your home directory.

Confguration files have to be updated changing /user/ to your login 
in keys which set paths. These keys are :
- 454_RUN_IDENTIFIER
- INPUT_FILE_FORWARD
- QUALITY_FILE or INPUT_FILE_REVERSE
- BARCODE_FILE
- SEQUENCE_FOLDER
- QUERY_SEQUENCES
- OTUDIST
- PHYLO_FOLDER
- CLADE_FILE
- REFERENCE_BASE

If you haven't installed PANAM in location said above, you also have to 
update the key REFERENCE_BASE, changing /home/user/panam/ to the real
PANAM's installation location.

Please pay attention to the value of the key FUZZNUC_PATH which has to
be the value returned by a command : which fuzznuc (if fuzznuc was
installed on your computer, otherwise install EMBOSS-6.5.7 which embed 
fuzznuc).
