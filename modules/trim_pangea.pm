#!/usr/bin/perl -w
#use strict;
#use warnings;

#############################
		# trim pangea  #
#############################


sub trim_pangea {
	## inputs
	my $MinSeqLength;
	my $MaxSeqLength;
	my $Nbase; 
	my $inputSeqNameF; 	 	
	my $NGS_id_Results;
	my $inputQualName;

	($inputSeqNameF, $MinSeqLength, $MaxSeqLength, $Nbase, $NGS_id_Results, $inputQualName) = @_ ;

	## outputs
	
	my $fasta_dir =  $NGS_id_Results."/fasta_dir/";

	if (-d $fasta_dir) {
		`rm -r $fasta_dir`;
		`mkdir $fasta_dir`
	}
	else {
		`mkdir $fasta_dir`
	}

	my $sampleAll = $fasta_dir."/seqAll_checkPrimers.fasta";
	my $tmp = $fasta_dir."/seqAll_tmp";

	## pangea_trim

	open(INPUTSEQ, $inputSeqNameF);
	open(INPUTQUAL, $inputQualName);
	open OUTPUT, ">$tmp" or die $!; 

	my $Rejected = -1;
	my @FinalTrim = ();
	my $end = 0;
	my $max = 0;
	my $sum = 0;
	my $start = 0;
	my $first = 0;
	my $lineNum = 0;

	
	while(my $lineSeq = <INPUTSEQ>)
	{
		my $lineQual = <INPUTQUAL>; $lineQual=~s/^\s+//;
		if($lineSeq =~ />/)
		{
			my $length = ($end + 1) - $start;
			if($length < $MinSeqLength)
			{
				$Rejected = 1;
			}
			if($Rejected == 0)
			{
				print OUTPUT "$header"; 
				print OUTPUT "\n";
				$countDown = 60;
				for($a = $start; $a <= $end; $a++)
				{
					$countDown--;
					if (defined $FinalTrim[$a]) {print OUTPUT "$FinalTrim[$a]";}
					if($countDown == 0)
					{
						print OUTPUT "\n";
						$countDown = 60;
					}
				}
				print OUTPUT "\n";
			}
			$Rejected = 0;
	
			#Trim off heading

	#############################
			if ($lineSeq=~ / /) {$space = index($lineSeq, " ") + 1}
			elsif ($lineSeq=~ /\n/) {$space = index($lineSeq, "\n") + 1}
	#############################	  

			$header = substr($lineSeq, 0, $space);
			chomp($header);
			@FinalTrim = ();
			$max = 0;
			$sum = 0;
			$first = 0;
			$lineNum = 0;
		}
		elsif($Rejected == 0)
		{
			#Store for use when trimming bases
			@bases = split(//, $lineSeq);
			$size  = scalar @bases;
			$size--;
			for($a = 0; $a < $size; $a++)
			{
				push(@FinalTrim, $bases[$a]);
			}
	
			chomp($lineQual);
			@line = split(/ /, $lineQual);
			$size = scalar @line;
			for($a = 0; $a < $size; $a++)
			{ 
				$sum += ($line[$a] - $MinQualScore);
				if($sum > $max)
				{
					$max = $sum;
					$end = $a + (60*$lineNum);
					$start = $first;
				}
				if($sum < 0)
				{
					$sum = 0;
					$first = $a + (60*$lineNum);
				}
			}
			$lineNum++;
		}
	}

	close INPUTSEQ;
	close INPUTQUAL;
	close OUTPUT;

	##############################

	open (TMP, $tmp);
	my %seq; my $a;
	while (<TMP>) {
		if ($_=~/^>(.*?)\n/) {
			$a = $1
		}
		else {
			chomp($_);
			$seq{$a}.= $_
		}
	}
	close TMP;
	
	open (NBASE,  ">$sampleAll");
	foreach my $k (keys %seq) {
		my $l = length ($seq{$k});
		if ($l >= $MinSeqLength) {
			if ($Nbase eq "no") {
				print NBASE ">$k\n$seq{$k}\n"
			}
			else {
				if ($seq{$k}=~ /N/i) {}
				else {
					print NBASE ">$k\n$seq{$k}\n"
				}
			}
		}
	}
	close NBASE;		

	`rm $tmp`;
	return ($fasta_dir);
}

1;
