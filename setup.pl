#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

################################""

my $usage = qq~
Usage: perl -w $0 -USEARCH_path <PATH>
       PATH:
                the full path of USEARCH 

~;

my ($usearch, $bioperl_home) = ();
GetOptions(	'USEARCH_path=s'=>\$usearch) || die $usage;

die $usage unless ($usearch);

unless (-e $usearch) { die "$usearch does not exist $!\n"; }

############################usearch
my @t = split(/\//, $usearch);
my $mm = $#t;

my $path;
my $uu = $t[$mm];
for (my $i=0; $i<$mm; $i++) {
	$path .=$t[$i]."/"
}

if (($uu=~ /usearch/ ) or ($uu=~ /uclust/)) {
	`chmod 755 $path$uu > file`;
}
		
die "\n\n\t=> Usearch, ".$path.$uu.", appears to be empty!\n\n" if (-z $path.$uu);
die "\n\n\t=> Usearch, ".$path.$uu.", is not an executable file! insure to make it executable before you launch setup.pl.\n\n" if (!(-z "file"));
`rm file`;

my $us_version;
my $version = `$path$uu --version`;
if ($version=~/\w.*?\sv(\d.*?)\s/) {
	$us_version = $1;
	chomp ($us_version); 
}

my $current_dir = `pwd`; 
chomp($current_dir);

my $cho;

if (!(defined $us_version)) { die "\nWrong USEARCH path!\nPlease enter the full path\n\n"};
	
print "\n\n\tYou have a license to use USEARCH version $us_version\n";
if (($us_version ne "3.0.617")) {
	print "\n\tYour version may not be suitable for the running of PANAM.\n\tFor review only, we provide version 3.0.617.\n\tPANAM will be set on this version.\n\n";
	$cho = $current_dir."/bin/uclust3.0.617";
	$us_version = "3.0.617";
	print "\n\tDo you want to continue? [y/n]\n\n\t";
	my $rep = <STDIN>; chomp ($rep);

	if (($rep ne "y") and ($rep ne "n")) {
		print "\n\n\tWrong answer\n\tPlease enter y or n\n\t";
		$rep = <STDIN>;
		chomp ($rep)
	}
			
	if ($rep eq "n") {die "\n\tInstalling PANAM interrupted\n\n"} 
}

else {
	$cho = $path.$uu;
}

################### check if pandaseq is installed

#my $panda = `which pandaseq` ;
#if ($panda eq "") {die "\n\tPandaseq is not found .... Installing PANAM interrupted\n\n"} 


################## setting panam.pl

my $fast= $current_dir."/bin/./FastTree";
my $align = $current_dir."/bin/./hmmer-2.3.2/src/hmmalign";
my $build = $current_dir."/bin/./hmmer-2.3.2/src/hmmbuild";
my $modules = $current_dir."/modules";
my $rscripts = $current_dir."/R";

open (F, "$current_dir/.panam.txt") ;
my $text1;
while (<F>) {
	$text1 .=$_
}
close F;

$text1 =~ s/USEARCH_PATH/$cho/g;
$text1 =~ s/US_VERSION/$us_version/g;
$text1 =~ s/FASTTREE_PATH/$fast/g;
$text1 =~ s/HMMALIGN_PATH/$align/g;
$text1 =~ s/HMMBUILD_PATH/$build/g;

open (FILE, ">$current_dir/panam.pl");
print FILE $text1;
close FILE;

###################" setting preprocess.pl
open (P, "$current_dir/.preprocess.txt") ;
my $text2;
while (<P>) {
	$text2 .= $_
}
close P;

$text2 =~ s/USEARCH_PATH/$cho/g;
$text2 =~ s/MODULES_PATH/$modules/g ;

open (F2, ">$current_dir/preprocess.pl");
print F2 $text2;
close F2;

###################" setting check_chimeras.pm
open (M, "$modules/.check_chimeras.txt") ;
my $text3;
while (<M>) {
	$text3 .= $_
}
close M;

$text3 =~ s/USEARCH_PATH/$cho/g;

open (F3, ">$modules/check_chimeras.pm") ;
print F3 $text3;
close F3;

###################" setting quality.pl

my $fuzz = $current_dir."/bin/EMBOSS-6.5.7/emboss/./fuzznuc";

open (Q, "$current_dir/.quality.txt") ;
my $text4;
while (<Q>) {
	$text4 .= $_
}
close Q;

$text4 =~ s/MODULES_PATH/$modules/g ;
$text4 =~ s/FUZZNUC_PATH/$fuzz/g ;

open (F4, ">$current_dir/quality.pl");
print F4 $text4;
close F4;

###################" setting postprocess.pl
open (PP, "$current_dir/.postprocess.txt") ;
my $text5;
while (<PP>) {
	$text5 .= $_
}
close PP;

$text5 =~ s/RSCRIPTS/$rscripts/g ;

open (F5, ">$current_dir/postprocess.pl");
print F5 $text5;
close F5;


###################" setting phylodiv.pl
open (PH, "$current_dir/.phylodiv.txt") ;
my $text6;
while (<PH>) {
	$text6 .= $_
}
close PH;

$text6 =~ s/RSCRIPTS/$rscripts/g ;
$text6 =~ s/MODULES_PATH/$modules/g ;
$text6 =~ s/TEXT-CSV_PATH/$current_dir\/bin\/Text-CSV_XS-1.04/g ;

open (F6, ">$current_dir/phylodiv.pl");
print F6 $text6;
close F6;

#######################################################""
print "\n\n\tCompiling Bioperl...\n\t******* For non experimented users, we recommand to use the default setting while installing bioperl *******\n";
print "\t******* If the archive Bioperl-1.5.2_102.zip already exists, you may chose to replace it or not [A/N] *******\n";
print "\t******* After just press enter when a choice is to be made *******\n\tContinue? [y/n]\n" ;


my $rep = <STDIN>; chomp ($rep);
if (($rep ne "y") and ($rep ne "n")) {
	print "\n\n\tWrong answer\n\tPlease enter y or n\n\t";
	$rep = <STDIN>;
	chomp ($rep)
}
if ($rep eq "n") {die "\n\tInstalling PANAM interrupted\n\n"}
else {
	compile_bioperl();
}

compile_fasttree();
unless (-x "FastTree") {die "\tFastTree has not been compiled $!\n" ; }

compile_hmmer();
unless ((-x "hmmer-2.3.2/src/hmmalign") and (-x "hmmer-2.3.2/src/hmmbuild")) {die "\tHMMER has not been compiled $!\n" ; }


unless (-e "bioperl-1.5.2_102/Bio/TreeIO.pm") { die "\tModule TreeIO does not exist $!\n"; }


install_text() ;

compile_emboss() ;
unless (-x "EMBOSS-6.5.7/emboss/./fuzznuc") {die "\tFuzznuc has not been found $!\n" ; }

sub compile_bioperl {
	chdir "bin";
	system ("unzip bioperl-1.5.2_102.zip");
	chdir "bioperl-1.5.2_102" ;

	system ("perl Build.PL --install_base ../"); 
	system ("./Build test");
	system ("./Build install");
	chdir ("../");
}


sub compile_fasttree {
	print "\nCompiling FastTree ...\n";
	system ("gcc FastTree.c -DNO_SSE -lm -O3 -finline-functions -funroll-loops -Wall -o FastTree");
}

sub compile_hmmer {
	print "\nCompiling hmmer ...\n";
	system ("tar -xvzf hmmer-2.3.2.tar.gz");
	chdir ("hmmer-2.3.2");
	system ("./configure");
	system ("make");
	chdir ("../");
}

sub install_text {
	print "\nInstalling text-csv module ...\n";
	system ("tar -xvzf Text-CSV_XS-1.04.tgz");
	chdir ("Text-CSV_XS-1.04");
	system ("perl Makefile.PL PREFIX=$current_dir/bin/Text-CSV_XS-1.04");
	system ("make test");
	system ("make install");
	chdir ("../");
}

sub compile_emboss {
	print "\nInstalling fuzznuc ...\n";
	system ("tar -xvzf EMBOSS-6.5.7.tar.gz");
	chdir ("EMBOSS-6.5.7");
	system ("./configure --prefix=$current_dir/bin/EMBOSS-6.5.7 --without-x");
	system ("make");
	chdir ("../");
}


###############################""

print "\n\n\tPANAM had been successfully installed.\n\n";
