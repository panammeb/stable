REQUIREMENTS
----------------------------------------------------------------------------------------------
Requirements: Linux OS. 

PANAM was tested on ubuntu system 12.4 and 14.04

The following softwares are required by the PANAM package. They need to be downloaded and installed separately from PANAM. 
1. Perl 5 or later (www.perl.org) 
2. USEARCH (www.drive5.com). Tested versions are v1.1.579q; v3.0.617; v4.0.38 and v5.0.150. 
Other versions may not be suitable to the running of PANAM. 
3. PANDASEQ 
4. gcc. 
5. make. 
6. R. 
7. R packages: Vegan; Phyloseq; Picante; Mass.

The following softwares are included in the PANAM package, they will be installed when setting PANAM up. (We do not guarantee the running of PANAM with other versions) 
1. FastTree-2.1.3 
2. HMMER-2.3.2 
3- Bioperl-1.5.2

GETTING PANAM
----------------------------------------------------------------------------------------------

Use git to clone stable repository :

	git clone https://github.com/panammeb/stable
	cd stable/

Or, download the zip archive :

	wget https://github.com/panammeb/stable/archive/master.zip
	unzip master.zip
	cd stable-master/

INSTALLING PANAM 
----------------------------------------------------------------------------------------------

Untar reference files :

	tar -xjvf Reference.tar.bz2
	rm Reference.tar.bz2

In aim to test panam, you can use the uclust binary included in the panam distribution.

Install panam with the setup.pl script :

	perl -w setup.pl -USEARCH_path /absolute/path/to/uclust/binary

Choose default choice for all questions.
Ignore FastTree compiling warnings if occurred.




