#!/usr/bin/perl -w
use strict;
use warnings;

#############################
		# pandaseq  #
#############################


sub merge_pairs_v2 {
	## inputs
	my $MinSeqLength;
	my $MaxSeqLength;
	my $Nbase; 
	my $MinOverlap;
	my $MismatchOverlap;
	my $inputSeqNameF; 	 	
	my $inputSeqNameR;
	my $NGS_id_Results;
	my $merge;
	my $score;

	($inputSeqNameF, $inputSeqNameR, $MinSeqLength, $MaxSeqLength, $Nbase, $MinOverlap, $MismatchOverlap, $NGS_id_Results, $merge, $score) = @_ ;


	## split fastq files
	#my $fastq_split_dir =  $NGS_id_Results."/fastq_split_dir/";

	#if (-d $fastq_split_dir) {
	#	`rm -r $fastq_split_dir`;
	#	`mkdir $fastq_split_dir`
	#}
	#else {
	#	`mkdir $fastq_split_dir`
	#}

	#`split -l 300000 -d  --suffix-length=6 $inputSeqNameF $fastq_split_dir/fastqF_`;
	#`split -l 300000 -d  --suffix-length=6 $inputSeqNameR $fastq_split_dir/fastqR_`;
	
	#my @listefic = <$fastq_split_dir*> ;

	#my %h;my %path;
	#foreach my $file (@listefic) { 
	#	if ($file =~ /(\/.*)fastqF_(.*)/ ) { 
	#	my $id = $2; 
	#		$h{$id}{'F'} = $file;
	#	}
	#	if ($file =~ /(\/.*)fastqR_(.*)/ ) { 
	#		my $id = $2;
	#		$h{$id}{'R'} = $file;
	#	}
	#}

	### merge fastq reads

	my $fasta_dir =  $NGS_id_Results."/fasta_dir/";
	my $log_dir = $NGS_id_Results."/pandaseq_log_dir/";

	if (-d $fasta_dir) {
		`rm -r $fasta_dir`;
		`mkdir $fasta_dir`
	}
	else {
		`mkdir $fasta_dir`
	}

	if (-d $log_dir) {
		`rm  -r $log_dir`;
		`mkdir $log_dir`
	}
	else {
		`mkdir $log_dir`
	}

	## pandaseq

	my $panda_dir = "/usr/bin/pandaseq" ;
	my $command;

	if ($merge eq "perfect_match") {
		#foreach my $k (keys %h) { 
			if ($Nbase eq "yes" ) {
				#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
				$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
			 
			}

			if ($Nbase eq "no" ) {
				#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
				$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -C completely_miss_the_point:".$MismatchOverlap." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";

			}

			system($command) ;
		#}
	}

	elsif ($merge eq "fix_sequences") {
			#foreach my $k (keys %h) { 
				if ($Nbase eq "yes" ) {
					#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
					$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -N -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";

				}

				if ($Nbase eq "no" ) {
					#$command = $panda_dir."./pandaseq -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";
					$command = $panda_dir." -f ".$inputSeqNameF." -r ".$inputSeqNameR." -l ".$MinSeqLength." -L ".$MaxSeqLength." -o ".$MinOverlap." -t ".$score." -g ".$log_dir."seqAll_log.txt -w ".$fasta_dir."seqAll_checkPrimers.fasta";

				}

				#print "$command\n";
				system($command) ;
			#}
		}

	my @fasta_file = <$fasta_dir*>;
	foreach my $fasta_file (@fasta_file) { 
		`sed -i s/:/-/g $fasta_file 2> /dev/null`;
		`sed -i s/ /-/g $fasta_file 2> /dev/null`;
	}

	#`rm -r $fastq_split_dir`;
	return ($fasta_dir, $log_dir);

}

1;
