##################################################
# universal_grass_peps pipeline v1.3             #
# Rowan Mitchell 2023                            #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################


package Blastp;

use strict;
use warnings;
use Pipeline;

my $HITSUBDIR = "hitseqs/";
my $GMDIR = "/home/data/tropical_forage_grasses/grass_core/final/EnsemblPlants/pep_gene_models/";
my $CT_LIM=1000; # min num lines accepted for blast output
my @KEY_TH_DETAILS = ("score",  "pc_id", "coverage", "query", "subject", "exclusive"); # output defaults
my $DELIM_TH_DETAILS = ";"; 
	

sub run_blastp {
# called 	my ($blastpF, $nql) = Blastp::run_blastp($sp, "Oryza_sativa");
	my ($qsp, $dbsp) = @_;
	my $fullqF = Pipeline::gmfasta_filename(  $qsp );
	my $fulldb = Pipeline::gmfasta_filename( $dbsp );
	my $blast_tabF = blastp_filename($qsp, $dbsp); 
	# used for rel52
#	`blastp -db $fulldb -query $fullqF -evalue 1.e-5 -max_target_seqs 50 -outfmt 7 -out $blast_tabF -num_threads 12`;
	#
	#+ used for rel55
	# for orphan families in Ensembl: "We use the version 2.2.28, with the parameters: -seg no -max_hsps_per_subject 1 -use_sw_tback -num_threads 1."
	# only using 1 hsp greatly simplifies processing. May be ok for current use
	if (-e $blast_tabF) {
		print "blast result file already exists: $blast_tabF Delete to run new.\n";
	} else {
		`blastp -db $fulldb -query $fullqF -evalue 1.e-5 -max_target_seqs 50 -seg no -max_hsps 1 -outfmt 7 -out $blast_tabF -num_threads 12`;
	}
	#-
	die "blast output file $blast_tabF not found\n" unless (-e $blast_tabF);
	open (CHK, $blast_tabF) || die "Cannot open $blast_tabF\n";
	my $ct=0;
	while (<CHK>) {$ct++ if ($_ =~ /Query:/)}
	close CHK;
#	die "blast output file $blast_tabF only has $ct Query: lines\n" unless ($ct > $CT_LIM);
	return ($blast_tabF, $ct);	
}

sub check_hitdir {
# called Blastp::check_hitdir();
	my $hitdir = $BLASTPDIR . $HITSUBDIR;
	my @fs = `ls $hitdir`;
	die "Directory $hitdir not empty.\nDelete content before running prog as it appends files\n" if (@fs>0);
	return;
}


sub passgroups_hits {
# get hits and find groups that pass
# called my ($passgrpids_ref, $sp_grp_hit_ref) = 
#					Blastp::passgroups_hits($grpids_ref, $sps_ref, $blastpths_ref, $ng_min);
	my ($grpids_ref, $sps_ref, $blastpths_ref, $ng_min) = @_;
	my @sps = @$sps_ref;
	my %blastpths = %$blastpths_ref;
	my @passgrpids; # 
	my %sp_grp_hit;
	foreach my $grpid (@$grpids_ref) {
		my $nhit=0;
		for (my $i=0; $i<@sps; $i++) {
			my ($nohit, $thtxtref) = parse_hitdetail( $blastpths{$grpid}[$i] );
			unless ($nohit) {
				$sp_grp_hit{ $sps[$i]}{$grpid} = $thtxtref->{query};
				$nhit++; 
			}
		}
		next if ($nhit < $ng_min ); 
		push @passgrpids, $grpid; 
	}
	return (\@passgrpids, \%sp_grp_hit);
}

sub blastp_filename {
	my ($qsp, $dbsp) = @_;
	return "$qsp\_$dbsp\.blastp.txt";
}

sub hitseq_filename {
# called	my $hitseqF = Blastp::hitseq_filename($grpid); 
	my $grpid = shift();
	my $hitseqF = $BLASTPDIR . $HITSUBDIR . $grpid . "_hitseq.fa";
	return $hitseqF;
}

sub unassigned_filename {
# called	my $faF = Blastp::unassigned_filename($sp);
	my $sp = shift();
	my $faF = $BLASTPDIR . $HITSUBDIR . $sp . "_unassigned_hits.fa";
	return $faF;
}
	
sub get_grpid {
# called my $grp_id = Blastp::get_grpid( \@prots );
	my $protsref = shift;
	my @prots = sort @$protsref; 
	my $grpid = $prots[0];
	my $n = @prots - 1;
	$grpid .= "_plus$n" if ($n > 0);
	return $grpid;
}	

	
sub get_blastp_hits {
# called 	my $hits_ref = Blastp::get_blastp_hits( $sp , $os_len_ref);
	# find best score hit for every query (Sq) subject ($s) combination
    # returned hash has for every Sq $s : {score => , pc_id=> %identity, type => tophsp/multihsp, hsp_details => []}
	# hsp is just row of top hsp or joined rows of multi hsps
#	my $sp_or_filename  = shift();
	my $blast_tabF  = shift();
	my $os_len_ref = shift();
	my $os_is = shift() || "s"; # is Os subject (s; default) or query (q)
	# my $blast_tabF;
	# if ($sp_or_filename =~ m/\//) { # arg is filename
		# $blast_tabF = $sp_or_filename;
	# } else {                       # arg is species
		# $blast_tabF= blastp_filename( $sp_or_filename); 
	# }
	my %os_len = %{$os_len_ref};
	my $oslen; 
	my %q_s = %{ parse_blast($blast_tabF)};
	my $TOL = 2; # tolerance of overlapping residues 
	# output hash:
	my %hits;
	foreach my $q (keys %q_s) {
		if ($os_is eq "q") {
			die "Os length not found for query $q\n" unless (exists $os_len {$q});
			$oslen = $os_len {$q};
		}
		foreach my $s (keys %{$q_s{ $q}} ) {
		# ----- get all HSPs
			my @hsps = @{ $q_s{$q}{$s}};
			@hsps = sort {$$b[9] <=> $$a[9]} @hsps; # sort by descending scores
			my $th_score = $hsps[0][9];
			# ------------following is all to deal wih multi hsps scoring better than top hsp
			# ----- find consistent multiple HSPs
			my @c_hsps; # consistent hsps
			push (@c_hsps, $hsps[0]);
			my %keep_hsp;
			for (my $i=1; $i<@hsps; $i++) {
				$keep_hsp{$i}=1;
			}
			while (scalar(keys %keep_hsp)>0) {
				foreach my $i (sort keys %keep_hsp) {
					my @a1 = @{ $hsps[$i]}; # check this one for Pipeline::overlaps with all in c_hsps
					foreach my $c_hsp (@c_hsps) {
						my @a2 = @{ $c_hsp};
						if (Pipeline::overlap( $a1[4], $a1[5], $a2[4], $a2[5])>$TOL || 
						Pipeline::overlap( $a1[6], $a1[7], $a2[6], $a2[7])>$TOL || 
						# check q and s order of hsps are same 
						(($a2[4] - $a1[4])*($a2[6] - $a1[6])) < 0) {
							delete $keep_hsp{$i};
							last;
						}
					}
					if (exists $keep_hsp{$i}) {
						push (@c_hsps, $hsps[$i]); 
						delete $keep_hsp{$i};
					}
				}	
			}
			#-- calculate coverage of c_hsps
			my ($j1, $j2); # fields of HSP corresponding to Os coordinates
			if ($os_is ne "q") {
				die "Os length not found for subject $s\n" unless (exists $os_len {$s});
				$oslen = $os_len {$s};
				$j1 = 6; $j2 =7;
			} else {
				$j1 = 4; $j2 =5;
			}
			my $score;
			my $pc_id;
			my $coverage;
			my $hsp_details;
			if (! defined $c_hsps[0][0]) {
				warn "no HSPs for q $q s $s \n";
				$score=0;
				$pc_id=0;
				$coverage=0;
				$hsp_details=[];
			} else {
				$score=$c_hsps[0][9];
				my $t_al=$c_hsps[0][1]; 
				my $t_wp=$c_hsps[0][0] * $c_hsps[0][1]; # to calculate pc_id and coverage for multiple HSPs
				my $toslen=($c_hsps[0][$j2]-$c_hsps[0][$j1]+1);
				for (my $i=1;$i<@c_hsps;$i++) {
					$score += $c_hsps[$i][9];
					$t_wp += ($c_hsps[$i][0] * $c_hsps[$i][1]);
					$t_al += $c_hsps[$i][1];
					$toslen += ($c_hsps[$i][$j2]-$c_hsps[$i][$j1]+1) # correct for any overlap with prev HSP
							   - Pipeline::overlap( $c_hsps[$i-1][$j1], $c_hsps[$i-1][$j2], $c_hsps[$i][$j1], $c_hsps[$i][$j2]);
				}
				$pc_id = ($t_wp / $t_al) * ($toslen / $oslen); # weighted %id of all hsps corr for unaligned len of Os
				$coverage = 100 * $toslen / $oslen;
				# ----- define hsp_details
				my @multihsp; # join fields of multihsps 
				my $nhspf = scalar @{ $c_hsps[0]}; # number of hsp fields
				for (my $j=0; $j<$nhspf; $j++) {
					my @cs;
					for (my $i=0; $i<@c_hsps; $i++) {
						if (defined $c_hsps[$i][$j]) {
							push @cs, $c_hsps[$i][$j];
						} else {
							push @cs, "";
						}
					}
					$multihsp[$j] = join(",", @cs);
				}
				$hsp_details = [ @multihsp ];
			}				
			#---------- put in output stucture
			$hits{ $q }{ $s } = {score=>$score,
			                     pc_id=>$pc_id,
								 coverage =>$coverage,
								 hsp_details=>$hsp_details};
		}
	}
	return \%hits;
}
	
	
sub parse_blast {	
# called 	my %q_s = %{ parse_blast($blast_tabF)};
	my $blastF = shift();
	die "File $blastF not found\n" unless (-e $blastF);
	open (IN, $blastF) || die "Cannot open input file $blastF\n";
	my %q_s; # store all blast results  
	while (<IN>) { 
		next if (/^#/);
		chomp;
	# Fields:  query id,  subject id,
	# @a  0 % identity, 1 alignment length, 2 mismatches, 3 gap opens, 
	# 4 q. start, 5 q. end, 6 s. start, 7 s. end, 8 evalue, 9 bit score
		my ($q, $s, @a) = split("\t",$_);
		# ----- save multiple HSPs
		my @hsps;
		@hsps = @{ $q_s{$q}{$s} } if (exists $q_s{$q}{$s});
		push (@hsps, [ @a ]);
		$q_s{$q}{$s} = [ @hsps];
	}
	close IN;
	return \%q_s;
}

sub get_subject_sets_exchits {
# called	my %grp_exchits = %{ Blast_hit::get_subject_sets_exchits ($hits_ref, \%pep_grp) };
	my ($ref_hits, $ref_s_sset) = @_;
	my %hits = %{ $ref_hits };
	my %s_sset = %{ $ref_s_sset };
	# to find all exclusive hits hit within a blast hit set for a list of subjects (e.g. group of similar prots)
	#+ first find tophit subjects for each query
	my %qth_s; # store subject which gives top score for this query
	foreach my $q (keys %hits) {
		foreach my $s (keys %{ $hits{ $q}} ) {
			die "? blastp subject $s not found in supplied groups\n" unless (exists $s_sset{ $s } );
			if (exists $qth_s{$q}) {
				my $ths = $qth_s{$q};
				$qth_s{$q} = $s if (_hit1_better_hit2( $hits{$q}{$s}, $hits{$q}{$ths}) );
			} else {
				$qth_s{$q} = $s;
			}
		}
	}
	#- 
	#+ find exc hits
	my %sset_hits; #output
	foreach my $q (keys %hits) {
		my $ths = $qth_s{$q};		
		foreach my $s (keys %{ $hits{ $q}} ) {
			my $sset_id = $s_sset{ $s };
			push @{ $sset_hits{$sset_id}{$s} }, $q if ($s_sset{ $ths } eq $sset_id ); # top hit subj is in same set
		}
	}
	#-
	#+ sort queries for each subject by descending score
	foreach my $sset_id (keys %sset_hits) {
		foreach my $s (keys %{ $sset_hits{ $sset_id }}) {
			my @qs = sort {$hits{$b}{$s}{score} <=> $hits{$a}{$s}{score}} @{ $sset_hits{ $sset_id }{$s} }; 
			$sset_hits{ $sset_id }{$s} = [@qs];
		}
	}
	return \%sset_hits;
}




sub get_subject_sets_tophits {
# called	$sp_osgrp_tophit{$sp} = Blast_hit::get_subject_sets_tophits ($hits_ref, \%prot_osgrp);
	my ($ref_hits, $ref_s_sset) = @_;
	my %hits = %{ $ref_hits };
	my %s_sset = %{ $ref_s_sset };
	my %sset_tophit ; # output
	# to find top hit within a blast hit set for a list of subjects (e.g. group of similar prots)
	my %qth_s; # store subject which gives top score for this query
	foreach my $q (keys %hits) {
		foreach my $s (keys %{ $hits{ $q}} ) {
			if (exists $qth_s{$q}) {
				my $ths = $qth_s{$q};
				$qth_s{$q} = $s if (_hit1_better_hit2( $hits{$q}{$s}, $hits{$q}{$ths}) );
			} else {
				$qth_s{$q} = $s;
			}
			die "? blastp subject $s not found in supplied groups\n" unless (exists $s_sset{ $s } );
			my $sset_id = $s_sset{ $s };
			next if (exists $sset_tophit{ $sset_id } && 
			         _hit1_better_hit2($sset_tophit{ $sset_id }, $hits{$q}{$s}) );
			$sset_tophit{ $sset_id }{query} = $q;
			$sset_tophit{ $sset_id }{subject} = $s;		
			foreach my $k (keys %{ $hits{$q}{$s}}) {
				$sset_tophit{ $sset_id }{$k} = $hits{$q}{$s}{$k};
			}
		}
	}
	# now define whether hit is exclusive or whether query matches better to subjects outside set
	foreach my $sset_id (keys %sset_tophit) {
		my $q = $sset_tophit{ $sset_id }{query};
		my $ths = $qth_s{$q};		
		if ($s_sset{ $ths } eq $sset_id ) {
			$sset_tophit{ $sset_id }{exclusive} = "yes";
		} else {
			my $qthscore = $hits{$q}{$ths}{score};
			$sset_tophit{ $sset_id }{exclusive} = "no matches better to $ths with $qthscore";
		}
	}
	return \%sset_tophit;
}

sub _hit1_better_hit2 {
	my ($h1, $h2) = @_;
	my $result = 0;
	if ($h1->{score} > $h2->{score}) {
		$result = 1;
	} elsif ($h1->{score} == $h2->{score}) {
		$result = 1 if ($h1->{pc_id} > $h2->{pc_id});
	}
	return $result;
}
	

sub write_hitdetail {
# called: 			my $cell = Blastp::write_hitdetail( $sp_osgrp_tophit{ $sp }{ $grpid } );			
	my %th = %{ shift() };
	my @txts;
	foreach my $k (@KEY_TH_DETAILS) {
		next unless (exists $th{ $k });
		if ($k eq "pc_id" || $k eq "coverage") {
			my $pctxt = sprintf ("%.2f%%", $th{ $k });
			push @txts, "$k=" . $pctxt;
		} else {
			push @txts, "$k=" . $th{ $k };
		}
	}
	return join($DELIM_TH_DETAILS, @txts);
}

sub read_tophits {
#called my ($grpids_ref, $grpid_osprots_ref, $sps_ref, $blastpths_ref) = Blastp::read_tophits($BLASTPF);
	my $thF = shift;
	my $inF = $BLASTPDIR . $thF;
	die "file $inF does not exist\n" unless (-e $inF);
	open (IN, $inF)|| die "Cannot open input file $inF\n";
	my $ct_row = <IN>;
	chomp $ct_row; # column title row
	my ($ct1, $ct2, @sp) = split ("\t", $ct_row);
	my @grpids;
	my %grpid_osprots;
	my %th;
	while (<IN>) { 
		chomp;
		my ($grpid, $protstxt, @cells) = split ("\t", $_); 
		$grpid_osprots{$grpid} = [split (" ", $protstxt)];
		push (@grpids, $grpid);
		$th { $grpid } = [ @cells ];
	}
	close IN;	
	return (\@grpids, \%grpid_osprots, \@sp, \%th);
}

sub parse_hitdetail {
# called: my ($nohit, $thtxtref) = Blast_hit::parse_hitdetail( $hitdets[ $iAt ] );
	my $cell = shift();
	my $nohit;
	my @a;
	$nohit = (! defined $cell || length( $cell) ==0);
	@a = split ($DELIM_TH_DETAILS, $cell) unless ($nohit);
	my %th_txt;
	foreach my $f (@a) {
		my ($k, $v) = split "=", $f;
		$th_txt{$k} = $v;
	}
	return ($nohit, \%th_txt);
}
1;
