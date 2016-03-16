REQUIREMENTS
----------------------------------------------------------------------------------------------
Requirements: Linux OS. 

PANAM was tested on ubuntu system 12.4 and 14.04

The following softwares are required by the PANAM package. They need to be downloaded and installed separately from PANAM. 

	1. Perl 5 or later (www.perl.org) 
	2. PANDASEQ 
	3. gcc 
	4. make 
	5. R 
	6. R packages: Vegan; Phyloseq; Picante; MASS

The following softwares are included in the PANAM package, they will be installed when setting PANAM up. (We do not guarantee the running of PANAM with other versions) 

	1. FastTree-2.1.3 
	2. HMMER-2.3.2 
	3. Bioperl-1.5.2
	4. UCLUST_v3.0.617
	5. PERL package: Text::CSV_XS
	6. fuzznuc included in EMBOSS-6.5.7

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

Install panam with the setup.pl script :

	perl -w setup.pl

Choose default choice for all questions.
Ignore FastTree compiling warnings if occurred.

Extract Reference archive :

	tar -xjvf Reference.tar.bz2


