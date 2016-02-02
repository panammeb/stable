#!/usr/bin/perl -w
use strict;
use warnings;

#############################
# Getting input parameters  #
#############################


sub parse_ini {

	my $option_file;
	($option_file) = @_ ;
	open(OFILE, "$option_file") || die "Cannot open input file : $option_file\n";

	######################################  Checks quality (Ns, min and max length, scores). Merges reads for illumina (merge,MinOverlap, MissmatchOverlap, score). Generates fasta files.

	my $NGS_id_Results;
	my $dataType;      	 	## 454 Roche/illumina MiSeq

	my $inputSeqNameF; 	 	## fasta/fastq_forward
	my $inputSeqNameR; 	 	## fastq_reverse ## if dataType= IlluminaMiSeq
	my $inputQualName;   	## if dataType= 454Roche 

	my $MinSeqLength;		## $lengthCut
	my $MaxSeqLength;		## $maxlength
	my $MinOverlap;			## if dataType= IlluminaMiSeq 
	
	my $Nbase;
	my $MinQualScore;	 	## if dataType= 454 only # $qualityCut

	my $merge;				## if dataType= IlluminaMiSeq
	my $score;				## if dataType= IlluminaMiSeq and merge= fix_sequence
	my $MismatchOverlap;	## if dataType= illuminaMiSeq and merge= perfect_match

	#####  Demultiplexing. Checks for primers and barcodes. Sorts sequence
	
	my $barIn;
	my $path_fuzznuc;
	my $forward;			## nom du primer F pour trimming profiles !!!
	my $reverse;			## nom du primer R pour trimming profiles !!!
	my $pmismatch ;
	my $match;


	###################################### Clustering and preprocessing

	my $seqFolder;		## Dossier des fichiers fasta à clusteriser fasta_dir
	my $Clst;
	my $nbr_seq_norm;	## seuil de normalisation
	my $pool;		### modif 05/03/2014
	

	###################### ajout 06/05/2014 ## variables pour phylodiv
	
	my $otuDist;
	my $clade;
	my $phyloFolder;
	my $refBase;

	###################### ajout 21/8/2014 ## variables pour CheckChimeras et CheckHomopol

	my $CheckChimeras;
	my $CheckHomopol;


	###################### ajout 21/10/2014 ## variable pour RemoveSingletons

	my $RemoveSingletons;

	############################################################################


	while (<OFILE>) {
		chomp;
		my $line = $_;

	# discard blank line
		if ($line =~ m/^\s*$/g) {
			next;
		}

	# discard comment line
		elsif($line =~ m/^\s*#/g) {
			next;
		}

	# get input options line
		else {

			if($line=~m/454_RUN_IDENTIFIER/gi) { 
				$line=~s/454_RUN_IDENTIFIER//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$NGS_id_Results = $line; 
				}
			}
	
			if($line=~m/DATA_TYPE/gi) { 
				$line=~s/DATA_TYPE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$dataType = $line; 
				}
			}

			if($line=~m/INPUT_FILE_FORWARD/gi) {
				$line=~s/INPUT_FILE_FORWARD//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$inputSeqNameF = $line;
				}
			}

			if($line=~m/INPUT_FILE_REVERSE/gi) {
				$line=~s/INPUT_FILE_REVERSE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$inputSeqNameR = $line;
				}
			}

			if($line=~m/QUALITY_FILE/gi) {
				$line=~s/QUALITY_FILE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$inputQualName = $line; 
				}
			}

			if($line=~m/LOWER_SEQUENCE_LENGTH_CUTOF/gi) {
				$line=~s/LOWER_SEQUENCE_LENGTH_CUTOF//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$MinSeqLength = $line; 
				}
			}

			if($line=~m/MAX_SEQUENCE_LENGTH/gi) {
				$line=~s/MAX_SEQUENCE_LENGTH//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$MaxSeqLength = $line; 
				}
			}

		
			if($line=~m/MIN_OVERLAP_LENGTH/gi) {
				$line=~s/MIN_OVERLAP_LENGTH//gi; 
				$line=~s/\t//gi; 
				$line=~s/^\s+//gi; 
				$line=~s/\s+$//gi; 
				if($line) {
					$MinOverlap = $line; 
				} 
			}
	
			if($line=~m/MISMATCH_OVERLAP/gi) {
				$line=~s/MISMATCH_OVERLAP//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi; 
				$MismatchOverlap = $line;
			}

			if($line=~m/AMBIGUOUS/gi) {
				$line=~s/AMBIGUOUS//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi; 
				if($line) {
					$Nbase = $line; 
				}
			}

			if($line=~m/AVERAGE_CUTOFF_QUALITY_VALUE/gi) {
				$line=~s/AVERAGE_CUTOFF_QUALITY_VALUE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$MinQualScore = $line; 
				}
			}

			if($line=~m/MERGE_SEQ/gi) {
				$line=~s/MERGE_SEQ//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$merge = $line; 
				}
			}	
		
			if($line=~m/SCORE_FIX/gi) {
				$line=~s/SCORE_FIX//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				$score = $line; 
			}	

			if($line=~m/BARCODE_FILE/gi) {
				$line=~s/BARCODE_FILE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$barIn = $line; 
				}
			}

			if ($line=~m/FUZZNUC_PATH/gi) {
				$line=~s/FUZZNUC_PATH//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi; 
				if($line) {
					$path_fuzznuc = $line; 
				}
			}

			if($line=~m/FORWARD_PRIMER_SEQUENCE/gi) {
				$line=~s/FORWARD_PRIMER_SEQUENCE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$forward = $line; 
				}
			}
		
			if($line=~m/REVERSE_PRIMER_SEQUENCE/gi) {
				$line=~s/REVERSE_PRIMER_SEQUENCE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$reverse = $line; 
				}
			}
		
			if($line=~m/PMISMATCH/gi) {
				$line=~s/PMISMATCH//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				$pmismatch = $line; 
			}

			if($line=~m/MATCH_PRIMER/gi) {
				$line=~s/MATCH_PRIMER//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$match = $line;  
				}
			}

			if($line=~m/SEQUENCE_FOLDER/gi) {
				$line=~s/SEQUENCE_FOLDER//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$seqFolder = $line; 
				}
				chomp ($seqFolder);
			}

			if($line=~m/CLUSTERING_CUTOFF/gi) {
			$line=~s/CLUSTERING_CUTOFF//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$Clst = $line; 
				}
			}

			if($line=~m/NBR_SEQ_NORM/gi) {
				$line=~s/NBR_SEQ_NORM//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if($line) {
					$nbr_seq_norm = $line;
				} 	
			}
			
			######### ajout 05/03/2014
			if($line=~m/POOLING/gi) {
				$line=~s/POOLING//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$pool = $line; 
				}
			}
			######### fin ajout 

			######### ajout 06/05/2014
			if ($line=~m/OTUDIST/gi) {
				$line=~s/OTUDIST//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$otuDist = $line; 
				}
			}

			if ($line=~m/CLADE_FILE/gi) {
				$line=~s/CLADE_FILE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$clade = $line; 
				}
			}

			if ($line=~m/PHYLO_FOLDER/gi) {
				$line=~s/PHYLO_FOLDER//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$phyloFolder = $line; 
				}
			}			
				
			if ($line=~m/REFERENCE_BASE/gi) {
				$line=~s/REFERENCE_BASE//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$refBase = $line; 
				}
			}

			##### ajout 21/8/2014
			if ($line=~m/CHECK_CHIMERAS/gi) {
				$line=~s/CHECK_CHIMERAS//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$CheckChimeras = $line; 
				}
			}	

			if ($line=~m/CHECK_HOMOPOLYMERS/gi) {
				$line=~s/CHECK_HOMOPOLYMERS//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$CheckHomopol = $line; 
				}
			}

			##### ajout 21/8/2014
			if ($line=~m/REMOVE_SINGL/gi) {
				$line=~s/REMOVE_SINGL//gi;
				$line=~s/\t//gi;
				$line=~s/^\s+//gi;
				$line=~s/\s+$//gi;
				if ($line) {
					$RemoveSingletons = $line; 
				}
			}		
		}
	}

	### ajout $pool ## 05/03/2014
	### ajout $otuDist, $clade, $phyloFolder, $refBase ## 06/05/2014
	### ajout $checkChimeras, $checkHomopol ## 21/8/2014
	### ajout $RemoveSingletons ## 21/10/2014

	return ($NGS_id_Results, $dataType, $inputSeqNameF, $inputQualName, $MinSeqLength, $MaxSeqLength, $MinOverlap, $MismatchOverlap, $Nbase, $MinQualScore, $inputSeqNameR, $merge, $score, $barIn, $path_fuzznuc, $forward, $reverse, $pmismatch, $match, $seqFolder, $Clst, $nbr_seq_norm, $pool, $otuDist, $clade, $phyloFolder, $refBase, $CheckChimeras, $CheckHomopol, $RemoveSingletons) 
}

1;
