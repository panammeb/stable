PANAM: Phylogenetic Analysis of Next generation AMplicons
-----------------------------------------------------------------------
TAIB N, BRONNER G, DEBROAS D.
Laboratoire Microorganismes : Génome et Environnement - UMR 6023 CNRS - FRANCE

-----------------------------------------------------------------------
Depict the microbial diversity using a phylogenetic approach

PANAM is a pipeline for depicting the microbial diversity within one or several samples. It implements a phylogenetic
approach for the annotation of 16S and 18S rRNA genes, and is optimized for the processing of large datasets
(Pyrosequencing and Illumina MiSeq runs). PANAM combines the publicly available tools: PANGEA, PANDASEQ,
USEARCH, HMMALIGN and FASTTREE with perl scripts dedicated to (i) cleaning raw sequences, (ii) trimming
profile alignments to fit the studied region and (iii) processing generated phylogenies to describe clades and
phylogenetic diversity indexes.
If you use this package, give this link as reference: https://github.com/panammeb/stable

-----------------------------------------------------------------------
License:
PANAM is free software: you may redistribute it and/or modify its under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2 of the License, or any later version.
PANAM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details
(http://www.gnu.org/licenses/).

-----------------------------------------------------------------------
Requirements:
Linux OS. PANAM was tested on ubuntu system 12.4 and 14.04
The following softwares are required by the PANAM package. They need to be downloaded and installed separately
from PANAM.
1. Perl 5 or later (www.perl.org)
2. USEARCH (www.drive5.com). Tested versions are v1.1.579q; v3.0.617; v4.0.38 and v5.0.150. Other versions may
not be suitable to the running of PANAM.
3. PANDASEQ
4. gcc.
5. make.
6. R.
7. R packages: Vegan; Phyloseq; Picante; Mass.

The following softwares are included in the PANAM package, they will be installed when setting PANAM up. (We do
not guarantee the running of PANAM with other versions)
1. FastTree-2.1.3
2. HMMER-2.3.2
3- Bioperl-1.5.2

-----------------------------------------------------------------------
