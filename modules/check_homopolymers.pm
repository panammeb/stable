#!/usr/bin/perl -w
#use strict;
#use warnings;


#############################
      # check homopolymers  #
#############################

sub check_homopolymers {
	my $sample_name;
	my $fasta_dir;
	($sample_name, $fasta_dir) = @_ ; 

	my $fasta_file = $fasta_dir."seqAll_".$sample_name.".fasta";
	my $minimum = 7;
	my $out_file = $fasta_dir."seqAll_".$sample_name."_homopolymers.fasta" ;
	
	if (!(-z $fasta_file)) {

	open(FASTA, "$fasta_file");
	open(OUT, ">$out_file") ;
	
	#####################################

	#---------------------------------------------------------------------------------------------------------------------------
	#The main event
	#---------------------------------------------------------------------------------------------------------------------------
	%headers = ();
	%records = ();
	$homopolymer_count = 9999;
	print "Poly Seq Length\tPolyA\tPolyT\tPolyG\tPolyC\tPolyN\tPolyX\n";
	@fasta = <FASTA>;
	$i = $minimum;
	while ($homopolymer_count > 0) {
	    $total_seq_length = 0;
	    $count = 0;
	    $counta = 0;
	    $countt = 0;
	    $countg = 0;
	    $countc = 0;
	    $countx = 0;
	    $countn = 0;
	    $string = "";
	    $header = "";
	    for $line (@fasta) {
		if ($line =~ /^>/) {
		    $count++;
		    $records{$header} = $string;
		    $homopolymer_count = &process_string($i);
		    $total_seq_length = $total_seq_length + length($string);
		    $header = $line;
		    chomp($header);
		    $header =~ s/\r//g;
		    $header =~ s/^>//g;
		    $string = "";
		}
		else {
		    chomp($line);
		    $string = $string . $line;
		}
	    }
	    $homopolymer_count = &process_string($i);
	    $total_seq_length = $total_seq_length + length($string);
	    $average_seq_length = $total_seq_length / $count;
	    $total_seq_length = 0;
	    print "$i\t$counta\t$countt\t$countg\t$countc\t$countn\t$countx\n";
	    $i++;
	}
	print "Number of Fasta Entries = $count\n";
	print "Average Sequence Length = $average_seq_length\n";
	close(FASTA);

	for $i (sort keys %headers) {
	    if ($headers{$i} =~ /\t/) {
		$acount = 0;
		$tcount = 0;
		$gcount = 0;
		$ccount = 0;
		$xcount = 0;
		$ncount = 0;
		@polycounts = split(/\t/, $headers{$i});
		for $j (0..$#polycounts) {
		    @parts = split(/=/, $polycounts[$j]);
		    if ($parts[0] eq "A" && $parts[1] > $acount) {
			$acount = $parts[1];
		    }
		    elsif ($parts[0] eq "T" && $parts[1] > $tcount) {
			$tcount = $parts[1];
		    }
		    elsif ($parts[0] eq "G" && $parts[1] > $gcount) {
			$gcount = $parts[1];
		    }
		    elsif ($parts[0] eq "C" && $parts[1] > $ccount) {
			$ccount = $parts[1];
		    }
		    elsif ($parts[0] eq "X" && $parts[1] > $xcount) {
			$xcount = $parts[1];
		    }
		    elsif ($parts[0] eq "N" && $parts[1] > $ncount) {
			$ncount = $parts[1];
		    }
		    else {
			print "ERROR!\n";
			&usage;
		    }
		}
		print OUT ">$i ";
		if ($acount > 0) {
		    print OUT "A=$acount";
		    if ($tcount > 0 || $gcount > 0 || $ccount > 0 || $xcount > 0 || $ncount > 0) {
			print OUT ", ";
		    }
		}
		if ($tcount > 0) {
		    print OUT "T=$tcount";
		    if ($gcount > 0 || $ccount > 0 || $xcount > 0 || $ncount > 0) {
			print OUT ", ";
		    }
		}
		if ($gcount > 0) {
		    print OUT "G=$gcount";
		    if ($ccount > 0 || $xcount > 0 || $ncount > 0) {
			print OUT ", ";
		    }
		}
		if ($ccount > 0) {
		    print OUT "C=$ccount";
		    if ($xcount > 0 || $ncount > 0) {
			print OUT ", ";
		    }
		}
		if ($xcount > 0) {
		    print OUT "X=$xcount";
		    if ($ncount > 0) {
			print OUT ", ";
		    }
		}
		if ($ncount > 0) {
		    print OUT "N=$ncount";
		}
		print OUT "\n";
		print OUT "$records{$i}\n";
	    }
	    else {
		print OUT ">$i $headers{$i}\n";
		print OUT "$records{$i}\n";
	    }
	}

	close(OUT);
	#--------------------------------------------------------------------------
	#Subroutines
	#--------------------------------------------------------------------------
	sub process_string {
	    ($myi) = @_;
	    $stringa = "A{$myi}";
	    $stringt = "T{$myi}";
	    $stringg = "G{$myi}";
	    $stringc = "C{$myi}";
	    $stringx = "X{$myi}";
	    $stringn = "N{$myi}";
	    if ($string =~ /$stringa/) { 
		$counta++;
		if ($headers{$header}) {
		    $headers{$header} = $headers{$header} . "\tA=" . $myi;
		}
		else {
		    $headers{$header} = "A=$myi";
		}
	    }
	    if ($string =~ /$stringt/) {
		$countt++;
		if ($headers{$header}) {
		    $headers{$header} = $headers{$header} . "\tT=" . $myi;
		}
		else {
		    $headers{$header} = "T=$myi";
		}
	    }
	    if ($string =~ /$stringg/) {
		$countg++;
		if ($headers{$header}) {
		    $headers{$header} = $headers{$header} . "\tG=" . $myi;
		}
		else {
		    $headers{$header} = "G=$myi";
		}
	    }
	    if ($string =~ /$stringc/) {
		$countc++;
		if ($headers{$header}) {
		    $headers{$header} = $headers{$header} . "\tC=" . $myi;
		}
		else {
		    $headers{$header} = "C=$myi";
		}
	    }
	    if ($string =~ /$stringx/) {
		$countx++;
		if ($headers{$header}) {
		    $headers{$header} = $headers{$header} . "\tX=" . $myi;
		}
		else {
		    $headers{$header} = "X=$myi";
		}
	    }
	    if ($string =~ /$stringn/) {
		$countn++;
		if ($headers{$header}) {
		    $headers{$header} = $headers{$header} . "\tN=" . $myi;
		}
		else {
		    $headers{$header} = "N=$myi";
		}
	    }
	    return $counta + $countt + $countg + $countc +$countx + $countn;
	}

	#####################################
	open (F, $fasta_dir."seqAll_".$sample_name."_homopolymers.fasta");
	my %h_homo; my $homo;
	while (my $l = <F>){
		chomp ($l);
		if ($l =~ m/^>(\S+)/){
			$homo = $1;
			$h_homo{$homo}++;
		}
	}
	close F;

	my %all_seq; my $seq_name;
	open (SEQ, $fasta_dir."seqAll_".$sample_name.".fasta"); 
	while (<SEQ>){
		chomp ($_);
		if ($_ =~ m/^>(.+)/){
			$seq_name = $1;
		}
		else {
			$all_seq{$seq_name} .= $_;
		}
	}
	close SEQ;

	open (R, ">".$fasta_dir."seqAll_".$sample_name."_homopolymersChecked.fasta");
	foreach my $kk (keys %all_seq) {	
		if (!(exists($h_homo{$kk}))){
			print R ">$kk\n$all_seq{$kk}\n";
		}
	}
	close R;

	`rm $fasta_dir"seqAll_"$sample_name"_homopolymers.fasta"`;
	#`mv $fasta_dir"seqAll_"$sample_name"_homopolymersChecked.fasta" $fasta_dir"/seqAll_"$sample_name".fasta"`;
	
}
}
1;
