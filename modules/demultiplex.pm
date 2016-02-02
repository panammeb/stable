#!/usr/bin/perl -w
#use strict;
#use warnings;


use check_chimeras;  ### ajout 21/8/2014
use check_homopolymers;   ### ajout 21/8/2014

#############################
	# demultiplex  #
#############################

sub demultiplex {

	my $NGS_id_Results;
	my $fasta_dir;
	my $fuzznuc;
	my $forward;
	my $reverse;
	my $pmismatch ;
	my $match;
	my $barIn;
	my $barcode; 
	my $MaxSeqLength;
	my $MinSeqLength;
		
	($NGS_id_Results, $fasta_dir, $fuzznuc, $forward,$reverse, $pmismatch, $match, $barIn, $barcode, $MaxSeqLength, $MinSeqLength , $CheckChimeras, $CheckHomopol) = @_ ; 
	#### ajout de $CheckChimeras 21/8/2014
	#### ajout de $CheckHomopol 21/8/2014

	my $tmp = $NGS_id_Results."/tmp";

	if (-d $tmp) {
		`rm -r $tmp`;
		`mkdir $tmp`;
	}
	else {`mkdir $tmp`}

	my $outputName = $tmp."/sample_all_oldIds_tmp";

	my @file_list= <$fasta_dir*> ;
	

	foreach my $file (@file_list) { 

		my $file_id;
		if ($file =~ /\/.*?fasta_dir\/(.*?).fasta/) {
			$file_id = $1
		}

		open (F, $file);		
		my %seq; my $a; 
		my %corresp; my $ii=1;
		my $id;
		my %ident;  ##???
	
		while(<F>){
			if ($_=~ /^>(.*?)\s+/) {
				$a =$1; 

				if (exists $ident{$a}) {
					 $a.="_".$ii;
				}
				else {
					$ident{$a} = 1   ###???
				}
				$id = "seq".$ii; 
				$corresp{$id} = $a;
				$ii++;	
			}
			else {
				$seq{$id}.=$_
			}
		}
		close F;

		my $corresp_file = $tmp."/".$file_id."_correspIds" ; 

		open (R, ">".$corresp_file) ;
		foreach my $k (keys %corresp) { 
			print R ">$k\n$seq{$k}\n"
		}
		close R;

		`$fuzznuc -sequence $corresp_file -pattern $forward -complement Y -pmismatch $pmismatch -outfile $tmp/$file_id"_outputF_tmp.fuzznuc" -rformat excel` ;
		`$fuzznuc -sequence $corresp_file -pattern $reverse -complement Y -pmismatch $pmismatch -outfile $tmp/$file_id"_outputR_tmp.fuzznuc" -rformat excel` ;
		

		#######

		my %fuzz_F; my %fuzz_R;

		my $i=0;
		my %posF; my %posR; 
	
		my %pos;
		my $startF; my $startR;	

		open (FUZF, $tmp."/".$file_id."_outputF_tmp.fuzznuc") ;
		open (RFF, ">".$tmp."/".$file_id."_outputF_Ids_tmp.fuzznuc") ;
		while (<FUZF>) {
			if ($_=~ /SeqName/) { print RFF "$_" }
			else {
				my @t = split (/\t/, $_); 
				$fuzz_F{$corresp{$t[0]}}++;
				print RFF "$corresp{$t[0]}\t$_";

				##### récupérer les positions de startF et startR pour les comparer, connaître le primer qui démarre le séquençage et chercher le tag en amont de cette amorce. 
				$pos{$corresp{$t[0]}}{'F'}{'start'} = $t[1]; 
				$pos{$corresp{$t[0]}}{'F'}{'end'} = $t[2];
 				if (exists $pos{$corresp{$t[0]}}{'F'}{'start'}) {$startF = $pos{$corresp{$t[0]}}{'F'}{'start'}}			
			}
		}	
		close FUZF;
		close RFF;

		open (FUZR, $tmp."/".$file_id."_outputR_tmp.fuzznuc") ;
		open (RFR, ">".$tmp."/".$file_id."_outputR_Ids_tmp.fuzznuc") ;
		my %h;
		while (<FUZR>) {
			if ($_=~ /SeqName/) { print RFR "$_"}
			else {
				my @t = split (/\t/, $_); 
				$fuzz_R{$corresp{$t[0]}}++;
				print RFR "$corresp{$t[0]}\t$_";
	
				###
				$pos{$corresp{$t[0]}}{'R'}{'start'} = $t[1]; ;
				$pos{$corresp{$t[0]}}{'R'}{'end'} = $t[2];
				if (exists $pos{$corresp{$t[0]}}{'R'}{'start'}) {$startR = $pos{$corresp{$t[0]}}{'R'}{'start'}}	
			}
		}	
		close FUZR;
		close RFR;

		if ($match eq "forward") {
			$i=0;
			open (FF, ">".$tmp."/".$file_id."_PrimersChecked_CheckBarcodes") || die "can not open file";
			foreach my $k (keys %seq) { 
				if (defined $fuzz_F{$corresp{$k}}) { 
					print FF ">$corresp{$k}\n$seq{$k}";
				}
			}
		}

		elsif ($match eq "reverse") {
			open (FF, ">".$tmp."/".$file_id."_PrimersChecked_CheckBarcodes") || die "can not open file";
			foreach my $k (keys %seq) {
				if (defined $fuzz_R{$corresp{$k}}) {
					print FF ">$corresp{$k}\n$seq{$k}";
				}
			}
		}

		elsif ($match eq "both") {
			open (FF, ">".$tmp."/".$file_id."_PrimersChecked_CheckBarcodes") || die "can not open file";
			foreach my $k (keys %seq) { 
				if ((exists $fuzz_F{$corresp{$k}}) and (exists $fuzz_R{$corresp{$k}})) { 
					print FF ">$corresp{$k}\n$seq{$k}";
				}
			}
		}

#		else { ## $match eq "none"
#			open (FF, ">".$tmp."/".$file_id."_PrimersChecked_CheckBarcodes") || die "can not open file";
#			foreach my $k (keys %seq) { 
#				print FF ">$corresp{$k}\n$seq{$k}";
#			}
#		}


		#############
			
		open (SI, $tmp."/".$file_id."_PrimersChecked_CheckBarcodes") || die "can not open file"; 
		my %seqq; my $aa;

		while (<SI>) {
		   if ($_=~ />(.*?)\s+/) { 
			  $aa = $1
		   }
			else {
			$seqq{$aa}.=$_
		   }
		}

		if ($barcode == 0 ) {
			open (SO, ">".$fasta_dir."/seqAll_oneSample.fasta"); 
			foreach my $kk (keys %seqq) {
				my $k = $corresp{$kk};	
				chomp ($seqq{$kk});			
				print SO ">oneSample_$kk\n$seqq{$kk}\n";
			}
			close SO;

			#### ajout 21/8/2014
			if (($CheckChimeras eq "yes") and ($CheckHomopol eq "no")) {
				my $sample_name = "oneSample";
				check_chimeras ($sample_name, $fasta_dir);

				`mv $fasta_dir"seqAll_"$sample_name"_chimerafree.fasta" $fasta_dir"SeqAll_"$sample_name".fasta"`;
			}
	
			if (($CheckChimeras eq "no") and ($CheckHomopol eq "yes")) {
				my $sample_name = "oneSample";
				check_homopolymers ($sample_name, $fasta_dir);

				`mv $fasta_dir"seqAll_"$sample_name"_homopolymersChecked.fasta" $fasta_dir"/seqAll_"$sample_name".fasta"`;
			}

			if (($CheckChimeras eq "yes") and ($CheckHomopol eq "yes")) {
				my $sample_name = "oneSample";
				check_chimeras ($sample_name, $fasta_dir);
			
				my $sample_name_chimeraChecked = $sample_name."_chimerafree";
				check_homopolymers ($sample_name_chimeraChecked, $fasta_dir);

				`rm $fasta_dir"seqAll_"$sample_name"_chimerafree.fasta"`;
				`mv $fasta_dir"seqAll_"$sample_name"_chimerafree_homopolymersChecked.fasta" $fasta_dir"seqAll_"$sample_name".fasta"`;			
			}
			#### fin ajout 21/8/2014
		}


		####
		my %sample_name;
		if ($barcode == 1) {
			my $l;
			my %barcode;

			open (B, $barIn) || die "can not open file"; 

			while (<B>) {
				my $line = $_;
				$line=~s/\s+$//gi;
				my @t = split (/\t/, $line); 
				chomp ($t[0]); chomp ($t[1]); 
				$barcode{$t[0]} = $t[1];
				$l = length($t[1]); 
			}
			close B;

			my %amont; my %res; my %posAbs;
			foreach my $k (keys %seqq){ 
				my $m; my $f;
				if ((defined $pos{$k}{'F'}{'start'}) and (defined $pos{$k}{'R'}{'end'})) {
					if ($startF < $startR) {	
						$m = $pos{$k}{'F'}{'start'} -1; 
						$f = $pos{$k}{'R'}{'end'}; 
					}
					else {
						$m = $pos{$k}{'R'}{'start'} -1;
						$f = $pos{$k}{'F'}{'end'}
					}
				}
				elsif ((defined $pos{$k}{'F'}{'start'}) and (!(defined $pos{$k}{'R'}{'end'}))) {
					$m = $pos{$k}{'F'}{'start'} -1;
				}
				#else {
				elsif ((!(defined $pos{$k}{'F'}{'start'})) and (defined $pos{$k}{'R'}{'end'})) {
					$m = $pos{$k}{'R'}{'start'} -1;
				}
	
				$posAbs{$k}{'start'} = $m;
				if (defined $f) {	
					$posAbs{$k}{'end'} = $f;	
#					print "m --- $m\nf --- $f\n";
				}
	
				my $mmoinsl = $m-$l;
				if ($mmoinsl >= 0) {
					$amont{$k}= substr($seqq{$k}, $mmoinsl, $l);  
					foreach my $bar (keys %barcode) { 
						if ($barcode{$bar} eq $amont{$k}) {  
							$res{$bar}{$k} = $seqq{$k}; 
						}
					}
				}
			}	

			my %clean_seq;


			foreach my $k (keys %res) { 
				foreach my $p (keys %{$res{$k}}) { 
					my $clean_seq;
					my $length;
					if (exists $posAbs{$p}{'end'}) {
						$length = $posAbs{$p}{'end'} - $posAbs{$p}{'start'}; 
					}
					else {
						$length = $MaxSeqLength
					}
		
					$clean_seq = substr ($seqq{$p}, $posAbs{$p}{'start'}, $length);
					my $ll = length ($clean_seq);
					if ($ll >= $MinSeqLength) {
						$clean_seq{$k}{$p}=$clean_seq;
					}
				}
			}
	
			#########################################"
	
			my $name ;
			if  ($file_id=~/(.*?)_checkPrimers/){
				$name = $1
			}
			foreach my $e (keys %clean_seq) { 
				open ($e, ">$fasta_dir/".$name."_".$e.".fasta") or die $!;
				foreach my $k (keys %{$clean_seq{$e}}) {
					my $l = length ($clean_seq{$e}{$k});
					chomp ($clean_seq{$e}{$k});
					print $e ">".$e."_"."$k\n$clean_seq{$e}{$k}\n"
				}
				close $e;
	
				##### ajout 21/8/2014
				if (($CheckChimeras eq "yes") and ($CheckHomopol eq "no")) {
					my $sample_name = $k;
					check_chimeras ($sample_name, $fasta_dir);
			
					`mv $fasta_dir"seqAll_"$sample_name"_chimerafree.fasta" $fasta_dir"SeqAll_"$sample_name".fasta"`;
				}
	
				if (($CheckChimeras eq "no") and ($CheckHomopol eq "yes")) {
					my $sample_name = $k;
					check_homopolymers ($sample_name, $fasta_dir);

					`mv $fasta_dir"seqAll_"$sample_name"_homopolymersChecked.fasta" $fasta_dir"/seqAll_"$sample_name".fasta"`;
				}

				if (($CheckChimeras eq "yes") and ($CheckHomopol eq "yes")) {
					my $sample_name = $k; 
					check_chimeras ($sample_name, $fasta_dir);

					my $sample_name_chimeraChecked = $sample_name."_chimerafree";
					check_homopolymers ($sample_name_chimeraChecked, $fasta_dir);

					`rm $fasta_dir"seqAll_"$sample_name"_chimerafree.fasta"`;
					`mv $fasta_dir"seqAll_"$sample_name"_chimerafree_homopolymersChecked.fasta" $fasta_dir"seqAll_"$sample_name".fasta"`;
				}
				##### fin ajout 21/8/2014

			}
		}
	}
	
	my @mv = <$fasta_dir*checkPrimers.fasta>;
	foreach my $e (@mv) {`mv $e $tmp`}

	return ($fasta_dir);

}

1;
		
