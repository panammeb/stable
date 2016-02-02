#!/usr/bin/perl -w
use strict;
use warnings;
#use fastx_barcode_splitter; 

use check_homopolymers;  ### ajout le 21/8/2014
use check_chimeras;   ### ajout le 21/8/2014

#############################
	# demultiplex_miseq #
#############################

sub demultiplex_miseq {

	my $NGS_id_Results;
	my $fasta_dir;
	my $fuzznuc;
	my $forward;
	my $reverse;
	my $pmismatch ;
	my $match;
	my $barIn;
	my $barcode; 
	my $MinSeqLength;
	my $CheckChimeras;
	my $CheckHomopol;
		
	($NGS_id_Results, $fasta_dir, $fuzznuc, $forward,$reverse, $pmismatch, $match, $barIn, $barcode, $MinSeqLength, $CheckChimeras, $CheckHomopol) = @_ ; 
	#### ajout de $MinSeqLength 10/6/2014 
	#### ajout de $CheckChimeras 21/8/2014
	#### ajout de $CheckHomopol 21/8/2014

	####
	my $tmp = $NGS_id_Results."/tmp";

	if (-d $tmp) {
		`rm -r $tmp`;
		`mkdir $tmp`;
	}
	else {`mkdir $tmp`}


	my @file_list= <$fasta_dir*> ; 

	my $fasta_split_dir =  $NGS_id_Results."/fasta_split_dir/";

	if (-d $fasta_split_dir) {
		`rm -r $fasta_split_dir`;
		`mkdir $fasta_split_dir`
	}
	else {
		`mkdir $fasta_split_dir`
	}

	####
	if ($barcode == 0 ) {

		open (SO, ">".$fasta_dir."seqAll_oneSample.fasta"); 

		## split fasta files
		
		foreach my $file (@file_list) {

			`split -l 300000 -d  --suffix-length=6 $file $fasta_split_dir"fasta_"`;		
		}
	

		my @split_list = <$fasta_split_dir*>;
		
		foreach my $split_file (@split_list) { 

			my $file_split_id;
			if ($split_file =~ /\/.*?fasta_split_dir\/(.*?fasta.*)/) {
				$file_split_id = $1; 
			}		

			open (FF, $split_file) or die "can not open file" ;
			my %seqq; my $aa;
			while (<FF>) {
				if ($_=~ />(.*?)\s+/) { 
			 		$aa = $1 ; 
		   		}
				else {
					$seqq{$aa}.=$_;
		 	 	}
			}
			close FF;
		

			`$fuzznuc -sequence $fasta_split_dir$file_split_id -pattern $forward -complement Y -pmismatch $pmismatch -outfile $tmp"/"$file_split_id"_outputF_tmp.fuzznuc" -rformat excel` ;
			`$fuzznuc -sequence $fasta_split_dir$file_split_id -pattern $reverse -complement Y -pmismatch $pmismatch -outfile $tmp"/"$file_split_id"_outputR_tmp.fuzznuc" -rformat excel` ;

			my %fuzz_F; my %fuzz_R;

			open (FUZF, $tmp."/".$file_split_id."_outputF_tmp.fuzznuc") ;
			while (<FUZF>) {
				if ($_=~ /SeqName/) { }
				else {
					my @t = split (/\t/, $_); 
					$fuzz_F{$t[0]}++;
				}
			}	
			close FUZF;

			open (FUZR, $tmp."/".$file_split_id."_outputR_tmp.fuzznuc") ;
			while (<FUZR>) {
				if ($_=~ /SeqName/) { }
				else {
					my @t = split (/\t/, $_); 
					$fuzz_R{$t[0]}++;
				}
			}	
			close FUZR;


			if ($match eq "forward") {
				open (FF, ">".$tmp."/".$file_split_id."_PrimersChecked") || die "can not open file";
				foreach my $k (keys %seqq) { 
					if (defined $fuzz_F{$k}) { 
						print FF ">$k\n$seqq{$k}";
						print SO ">oneSample_$k\n$seqq{$k}\n"; 
					}
				}
			}

			elsif ($match eq "reverse") {
				open (FF, ">".$tmp."/".$file_split_id."_PrimersChecked") || die "can not open file";
				foreach my $k (keys %seqq) {
					if (defined $fuzz_R{$k}) {
						print FF ">$k\n$seqq{$k}";
						print SO ">oneSample_$k\n$seqq{$k}\n";
					}
				}
			}

			elsif ($match eq "both") {
				open (FF, ">".$tmp."/".$file_split_id."_PrimersChecked") || die "can not open file";
				foreach my $k (keys %seqq) { 
					if ((exists $fuzz_F{$k}) and (exists $fuzz_R{$k})) { 
						print FF ">$k\n$seqq{$k}";
						print SO ">oneSample_$k\n$seqq{$k}\n";
					}
				}
			}

		}
		close SO;

		#### ajout de $CheckChimeras et CheckHomopol 21/8/2014

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
		############
	}

	#####
	if ($barcode == 1 ) {

		my $tag_tmp = $NGS_id_Results."/tag_tmp";

		if (-d $tag_tmp) {
			`rm -r $tag_tmp`;
			`mkdir $tag_tmp`;
		}
		else {`mkdir $tag_tmp`}

		my %dualbar;

		open (BAR, $barIn) || die "can not open file"; 
		while (<BAR>) {
			chomp($_);
			my @t = split (/\t/, $_);
			if ((exists $t[1]) and ($t[1] ne "")) {
				$dualbar{$t[0]}{'F'} = $t[1];		
			}
			if ((exists $t[2]) and ($t[2] ne "")) {
				$dualbar{$t[0]}{'R'} = $t[2]
			}
		}
		close BAR;

	
		#############################" modif 02/7/2014
		##################################################""

		## split fasta files
		
		foreach my $file (@file_list) {
			`split -l 300000 -d  --suffix-length=6 $file $fasta_split_dir"fasta_"`;		
		}
	

		my @split_list = <$fasta_split_dir*>;
			
		my %fuzz_F; my %rev_seq; my %seqq;		

		foreach my $split_file (@split_list) { 

			my $file_split_id;
			if ($split_file =~ /\/.*?fasta_split_dir\/(.*?fasta.*)/) {
				$file_split_id = $1; 
			}		

			open (FF, $split_file) or die "can not open file" ;
			my $aa;
			while (<FF>) {
				if ($_=~ />(.*?)\s+/) { 
			 		$aa = $1 ; 
					chomp($aa);
		   		}
				else {
					$seqq{$aa}.=$_;
		 	 	}
			}
			close FF;

			`$fuzznuc -sequence $fasta_split_dir$file_split_id -pattern $forward -complement Y -pmismatch $pmismatch -outfile $tmp"/"$file_split_id"_outputF_tmp.fuzznuc" -rformat excel` ;

			open (FUZF, $tmp."/".$file_split_id."_outputF_tmp.fuzznuc") ;
			while (<FUZF>) {
				if ($_=~ /SeqName/) { }
				else {
					my @t = split (/\t/, $_); 
					if ($t[4] eq "-") {
						$rev_seq{$t[0]} = reverse($seqq{$t[0]});
						$rev_seq{$t[0]} =~ tr/ACGTRMBD/TGCAYKVH/;
						chomp($rev_seq{$t[0]});					
						my $len = length($seqq{$t[0]});
						$fuzz_F{$t[0]}{'start'} = $len - $t[1] + 1;
						$fuzz_F{$t[0]}{'end'} = $fuzz_F{$t[0]}{'start'} + $t[3]; 
					}
					else {
						$fuzz_F{$t[0]}{'start'} = $t[1];
					}					
				}
			}	
			close FUZF;	
		}

		###########
		my %strpl;  ####### modif 25/06/2015

		open (STRPL, ">".$fasta_dir."seqAll_checkPrimers_strandplus");
		#open (STRPL, ">".$fasta_dir."seqAll_checkPrimers_strandplus_tmp");
		foreach my $aq (keys %seqq) {
			chomp($aq);
			if (exists $rev_seq{$aq}) {
				chomp($rev_seq{$aq});
				print STRPL ">$aq\n$rev_seq{$aq}\n";
				$strpl{$aq} = $rev_seq{$aq} ;  ##### modif 25/06/2015
		
			}
			#else {   #### modif 01/06/2015
			#	print STRPL ">$aq\n$seqq{$aq}\n";
			else {
				if ((exists $fuzz_F{$aq}{'start'}) and (exists $fuzz_F{$aq}{'start'} ne "")) {
					print STRPL ">$aq\n$seqq{$aq}\n";
					$strpl{$aq} = $seqq{$aq} ;  ##### modif 25/06/2015
				}
			}
		}
		close STRPL; 
		
		#############################
   
		#my $tmp_file = $fasta_dir."seqAll_checkPrimers_strandplus_tmp";
		#my $new_file = $fasta_dir."seqAll_checkPrimers_strandplus";
		#
		#open NOUVEAU, '>', $new_file;
		#open ANCIEN, '<', $tmp_file;
		#
		#while (<ANCIEN>) { # ligne stockée automatiquement dans $_;
		#	next if /^\s*$/;
		#	print NOUVEAU;
		#}
		#
		#close ANCIEN;
		#close NOUVEAU;
		#
		#`rm $tmp_file`;	
		
		################################

		foreach my $k (keys %dualbar) { print "k --- $k\n";
			my %fuzz_R;
			my %keep;

			my $file_to_split = $tmp."/".$k ;
			open (FTS, ">".$file_to_split) ;
			#foreach my $us (keys %seqq) {
			foreach my $us (keys %strpl) { ##### modif 25/06/2015
				#chomp($us); chomp($seqq{$us});
				chomp($us); chomp($strpl{$us});  ##### modif 25/06/2015
				#my $first30 = substr $seqq{$us}, 0, 30; 
				my $first30 = substr $strpl{$us}, 0, 30; 			
				if ($first30 =~ /$dualbar{$k}{'F'}/) {
					if (exists $dualbar{$k}{'R'}) {
						#my $last30 = substr $seqq{$us}, - 30;
						my $last30 = substr $strpl{$us}, - 30; ##### modif 25/06/2015
						my $rev = reverse($dualbar{$k}{'R'}); 
						$rev =~tr/ACGT/TGCA/; 
						if ($last30 =~ /$rev/) {
							#print FTS ">$us\n$seqq{$us}\n";
							print FTS ">$us\n$strpl{$us}\n";  ##### modif 25/06/2015
							$keep{$us}=1;
						}
					}
					else {
						#print FTS ">$us\n$seqq{$us}\n";
						print FTS ">$us\n$strpl{$us}\n";  ##### modif 25/06/2015
						$keep{$us}=1;
					}
				}
			}
		
			close FTS;	


############### spliter le résultat avant de vérifier amorces et couper les barcodes !!!
############" récupérer intersection barcodeF et barcodeR si les deux existent !!
				
			`split -l 300000 -d  --suffix-length=6 $file_to_split $fasta_split_dir$k".fasta_"`;	

			my @split_list = <$fasta_split_dir*$k.fasta*>;
		
			foreach my $split_file (@split_list) { 
				my $file_split_id;
				if ($split_file =~ /\/.*?fasta_split_dir\/(.*?fasta.*)/) {
					$file_split_id = $1; 
				}		

				#open (FF, $split_file) or die "can not open file" ;
				#my %seqq; my $aa;
				#while (<FF>) {
				#	if ($_=~ />(.*?)\s+/) { 
		 		#		$aa = $1 ; 
	   			#	}
				#	else {
				#		$seqq{$aa}.=$_;
	 	 		#	}
				#}
				#close FF;

				`$fuzznuc -sequence $fasta_split_dir$file_split_id -pattern $reverse -complement Y -pmismatch $pmismatch -outfile $tmp"/"$file_split_id"_outputR_tmp.fuzznuc" -rformat excel` ;
				open (FUZR, $tmp."/".$file_split_id."_outputR_tmp.fuzznuc") ;
				while (<FUZR>) {
					if ($_=~ /SeqName/) { }
					else {
						my @t = split (/\t/, $_); 
						#$fuzz_R{$t[0]}++;
						$fuzz_R{$t[0]}{'start'} = $t[1]; 
						$fuzz_R{$t[0]}{'end'} = $t[2];
					}
				}	
				close FUZR;

			}
		
			open (OS, ">".$fasta_dir."seqAll_".$k.".fasta");
			#foreach my $kk (keys %keep) {
			foreach my $kk (keys %fuzz_R) { 
				if ((exists $fuzz_F{$kk}{'start'}) and (exists $fuzz_R{$kk}{'end'})) { 
					my $m = $fuzz_F{$kk}{'start'} -1;
					my $length = $fuzz_R{$kk}{'end'} - $m ;
					my $seq = substr ($seqq{$kk}, $m, $length);
					######## modif 10/6/2014
					my $len = length $seq ;
					if ($len >= $MinSeqLength) {
						print OS ">".$k."_".$kk."\n$seq\n"
					}
					######### fin modif
				}
			}
			close OS;
		
			############ ajout $CkechChimeras et $CheckHomopol 21/8/2014
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
			############	
	
		}

#		##################################
#
#		my @file_list_strandplus =  <$fasta_dir*strandplus> ; 
#
#		foreach my $k (keys %dualbar) { print "k -- $k\n";
#	
#			open (DUF, ">".$tag_tmp."/".$k."_F") ;
#			print DUF "$k\t$dualbar{$k}{'F'}\n";
#			close DUF;
#
#			if (exists $dualbar{$k}{'R'}) {
#				open (DUR, ">".$tag_tmp."/".$k."_R") ;
#				print DUR "$k\t$dualbar{$k}{'R'}\n";
#				close DUR;
#			}
#
#			open (OS, ">".$fasta_dir."seqAll_".$k.".fasta");
#			my $file_to_split;
#
#			foreach my $file (@file_list_strandplus) { print "file --- $file\n";
#
#				fastx_split({"inputfile" => "$file", "bcfile" => "$tag_tmp/$k"."_F", "bol" => 1, "exact" => 1, "prefix" => "$tmp/", "suffix" => "_F" });
#				$file_to_split = $tmp."/".$k."_F";
#				
#				if (exists $dualbar{$k}{'R'}) {
#					fastx_split({"inputfile" => "$tmp/$k"."_F", "bcfile" => "$tag_tmp/$k"."_R", "eol" => 1, "exact" => 1, "prefix" => "$tmp/", "suffix" => "_R" });
#					$file_to_split = $tmp."/".$k."_R";	
#				}
#	
############### spliter le résultat avant de vérifier amorces et couper les barcodes !!!
############" récupérer intersection barcodeF et barcodeR si les deux existent !!
#				
#				`split -l 300000 -d  --suffix-length=6 $file_to_split $fasta_split_dir$k".fasta_"`;	
#
#				my @split_list = <$fasta_split_dir*$k.fasta*>;
#		
#				foreach my $split_file (@split_list) { 
#
#					my $file_split_id;
#					if ($split_file =~ /\/.*?fasta_split_dir\/(.*?fasta.*)/) {
#						$file_split_id = $1; 
#					}		
#
#					open (FF, $split_file) or die "can not open file" ;
#					my %seqq; my $aa;
#					while (<FF>) {
#						if ($_=~ />(.*?)\s+/) { 
#			 				$aa = $1 ; 
#		   				}
#						else {
#							$seqq{$aa}.=$_;
#		 	 			}
#					}
#					close FF;
#
#					`$fuzznuc -sequence $fasta_split_dir$file_split_id -pattern $reverse -complement Y -pmismatch $pmismatch -outfile $tmp"/"$file_split_id"_outputR_tmp.fuzznuc" -rformat excel` ;
#
#					open (FUZR, $tmp."/".$file_split_id."_outputR_tmp.fuzznuc") ;
#					while (<FUZR>) {
#						if ($_=~ /SeqName/) { }
#						else {
#							my @t = split (/\t/, $_); 
#							#$fuzz_R{$t[0]}++;
#							$fuzz_R{$t[0]}{'start'} = $t[1]; 
#							$fuzz_R{$t[0]}{'end'} = $t[2];
#						}
#					}	
#					close FUZR;
#
#					foreach my $kk (keys %seqq) { 
#						if ((exists $fuzz_F{$kk}{'start'}) and (exists $fuzz_R{$kk}{'end'})) { 
#							my $m = $fuzz_F{$kk}{'start'} -1;
#							my $length = $fuzz_R{$kk}{'end'} - $m ;
#							my $seq = substr ($seqq{$kk}, $m, $length);
#							######## modif 10/6/2014
#							my $len = length $seq ;
#							if ($len >= $MinSeqLength) {
#								print OS ">".$k."_".$kk."\n$seq\n"
#							}
#							######### fin modif
#						}
#					}
#				}
#				close OS;
#			}
#		}


		############################################## fin modif 2/7/2014
		###########################################		
#		
#		foreach my $k (keys %dualbar) { print "k -- $k\n";
#	
#			open (DUF, ">".$tag_tmp."/".$k."_F") ;
#			print DUF "$k\t$dualbar{$k}{'F'}\n";
#			close DUF;
#
#			if (exists $dualbar{$k}{'R'}) {
#				open (DUR, ">".$tag_tmp."/".$k."_R") ;
#				print DUR "$k\t$dualbar{$k}{'R'}\n";
#				close DUR;
#			}
#
#			open (OS, ">".$fasta_dir."seqAll_".$k.".fasta");
#			my $file_to_split;
#
#			foreach my $file (@file_list) { 
#
#				#`cat $file | perl /home/ntaib/panam_v4/fastx_barcode_splitter.pl --bcfile $tag_tmp/$k"_F" --bol --exact --prefix $tmp/ --suffix "_F" ` ;
#				fastx_split({"inputfile" => "$file", "bcfile" => "$tag_tmp/$k"."_F", "bol" => 1, "exact" => 1, "prefix" => "$tmp/", "suffix" => "_F" });
#				$file_to_split = $tmp."/".$k."_F";
#				
#				if (exists $dualbar{$k}{'R'}) {
#					#`cat $tmp/$k"_F" | perl /home/ntaib/panam_v4/fastx_barcode_splitter.pl --bcfile $tag_tmp/$k"_R" --eol --exact --prefix $tmp/ --suffix "_R" ` ;
#				fastx_split({"inputfile" => "$tmp/$k"."_F", "bcfile" => "$tag_tmp/$k"."_R", "eol" => 1, "exact" => 1, "prefix" => "$tmp/", "suffix" => "_R" });
#					$file_to_split = $tmp."/".$k."_R";	
#				}
#			
############### spliter le résultat avant de vérifier amorces et couper les barcodes !!!
############" récupérer intersection barcodeF et barcodeR si les deux existent !!
#				
#				`split -l 300000 -d  --suffix-length=6 $file_to_split $fasta_split_dir$k".fasta_"`;	
#
#				my @split_list = <$fasta_split_dir*$k.fasta*>;
#		
#				foreach my $split_file (@split_list) { 
#
#					my $file_split_id;
#					if ($split_file =~ /\/.*?fasta_split_dir\/(.*?fasta.*)/) {
#						$file_split_id = $1; 
#					}		
#
#					open (FF, $split_file) or die "can not open file" ;
#					my %seqq; my $aa;
#					while (<FF>) {
#						if ($_=~ />(.*?)\s+/) { 
#			 				$aa = $1 ; 
#		   				}
#						else {
#							$seqq{$aa}.=$_;
#		 	 			}
#					}
#					close FF;
#
#					`$fuzznuc -sequence $fasta_split_dir$file_split_id -pattern $forward -complement Y -pmismatch $pmismatch -outfile $tmp"/"$file_split_id"_outputF_tmp.fuzznuc" -rformat excel` ;
#					`$fuzznuc -sequence $fasta_split_dir$file_split_id -pattern $reverse -complement Y -pmismatch $pmismatch -outfile $tmp"/"$file_split_id"_outputR_tmp.fuzznuc" -rformat excel` ;
#		
#					my %fuzz_F; my %fuzz_R;
#
#					open (FUZF, $tmp."/".$file_split_id."_outputF_tmp.fuzznuc") ;
#					while (<FUZF>) {
#						if ($_=~ /SeqName/) { }
#						else {
#							my @t = split (/\t/, $_); 
#							#$fuzz_F{$t[0]}++;
#							$fuzz_F{$t[0]}{'start'} = $t[1]; 
#							$fuzz_F{$t[0]}{'end'} = $t[2];
#						}
#					}	
#					close FUZF;
#
#					open (FUZR, $tmp."/".$file_split_id."_outputR_tmp.fuzznuc") ;
#					while (<FUZR>) {
#						if ($_=~ /SeqName/) { }
#						else {
#							my @t = split (/\t/, $_); 
#							#$fuzz_R{$t[0]}++;
#							$fuzz_R{$t[0]}{'start'} = $t[1]; 
#							$fuzz_R{$t[0]}{'end'} = $t[2];
#						}
#					}	
#					close FUZR;
#					
#					foreach my $kk (keys %seqq) { 
#						if ((exists $fuzz_F{$kk}{'start'}) and (exists $fuzz_R{$kk}{'end'})) { 
#							my $m = $fuzz_F{$kk}{'start'} -1;
#							my $length = $fuzz_R{$kk}{'end'} - $m ;
#							my $seq = substr ($seqq{$kk}, $m, $length);
#							######## modif 10/6/2014
#							my $len = length $seq ;
#							if ($len >= $MinSeqLength) {
#								print OS ">".$k."_".$kk."\n$seq\n"
#							}
#							######### fin modif
#						}
#					}
#				}
#			}close OS;
#		}
########################################
	}
}

1;

























