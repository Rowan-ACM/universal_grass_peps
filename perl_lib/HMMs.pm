##################################################
# universal_grass_peps pipeline v1.3             #
# Rowan Mitchell 2023                            #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################


package HMMs;

use strict;
use warnings;
use Pipeline;


# column titles of hmmscan output files
my %COLTITLES = (hmmscan=>["target","target_acc","query","query_acc",
					 "Eval_full","score_full","bias_full",
					 "Eval_1dom","score_1dom","bias_1dom",
					 "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
					 "desc_target"], 
				hmmscandom=>["target","target_acc","tlen","query","query_acc","qlen",                   
					 "Eval_full","score_full","bias_full",
					 "domnum","of_doms","c-Evalue","i-Evalue","score_dom","bias_dom",
					 "from_hmm_coord","to_hmm_coord","from_ali_coord","to_ali_coord","from_env_coord","to_env_coord",
					 "acc","desc_target"] );

# column titles of hmmscan hit files with output for each query
my @HMMSCANHITSCOLS = ("query", "grpid", "hit_score", "maxposs_score", "hit_relscore", "lowestmem_relscore", "diff_lowestmem_hit_relscore");

		# hmmer files to keep by default: Full set is defined in Pipeline:
	# my %grpFs = (fasta=>$faF, mla=>$mlaF, hmm=>$hmmF, sto=>$stoF, hmmscan=>$hmmscanF, hmmscandom=>$hmmscandomF, 
	# consensus=>$confaF, conhmmscan=>$conhmmscanF, conhmmscandom=>$conhmmscandomF,
	# alt=>$altfaF, althmmscan=>$althmmscanF, althmmscandom=>$althmmscandomF);
my @KEEPDEF = ("fasta", "mla", "hmm", "sto", "consensus", "alt");

	
sub align_group {
# called	my $hmmFsR = HMMs::align_group($grpid, $hitseqF); 
# or	my $hmmFsR = align_group($grpid); 
	my $grpid = shift();
	my $useHMM0 = (my $HMM0F = shift());
	my $hmmFsR = Pipeline::group_fnames($grpid);
	my $faF = $hmmFsR->{fasta};
	`cp $HMM0F $faF`if ($useHMM0);
	die "make_groupdb: file $faF not there\n" unless (-e $faF);
	# step1 muscle
	my $mlaF = $hmmFsR->{mla};
	`muscle -in $faF -out $mlaF -quiet`;
	die "No muscle output file $mlaF created\n" unless (-e $mlaF);
	# step2 hmmbuild
	my $hmmF = $hmmFsR->{hmm};
	my $stoF = $hmmFsR->{sto};
	`hmmbuild --amino --fragthresh 0 -O $stoF $hmmF $mlaF`;
	die "No hmmbuild output file $stoF created\n" unless (-e $stoF);
	return $hmmFsR;
}

sub make_groupdb {
# called	$grpdetails_ref = HMMs::make_groupdb($grpid);
	my $hmmFsR = align_group( @_);
	my $hmmF = $hmmFsR->{hmm};
	# create db with hmmpress
	`hmmpress -f $hmmF`;
	# scan against self members with hmmscan
	my $faF = $hmmFsR->{fasta};
	my $hmmscanF = $hmmFsR->{hmmscan};
	_run_hmmscan($hmmFsR->{hmmscan}, $hmmFsR->{hmmscandom}, $hmmF, $faF);
	my %hits = %{ get_hmm_scores($hmmFsR) };
	#+ ------ Get max possible score from hmmscan of consensus
	my $confaF = $hmmFsR->{consensus};
	`hmmemit -c -o $confaF $hmmF`;
	my $conhmmscanF = $hmmFsR->{conhmmscan};
	my $conhmmscandomF = $hmmFsR->{conhmmscandom};
	_run_hmmscan($conhmmscanF, $conhmmscandomF, $hmmF, $confaF);
	my %conhmmFs = (hmmscan=>$conhmmscanF, hmmscandom=>$conhmmscandomF);
	my %conhits = %{ get_hmm_scores( \%conhmmFs) };
	my @qs = keys %conhits; # queries are ids for seqs in fasta file
	die "n queries found in consensus self hit file $conhmmscanF not = 1\n" unless (@qs == 1);
	my @ts = keys %{ $conhits{$qs[0]}};
	die "n targets found in consensus self hit file $conhmmscanF not = 1\n" unless (@ts == 1);
	my $maxposssc = $conhits{$qs[0]}{$ts[0]}{"best_score"};
	#- 
	#+ -----  get hmm relative self scores
	# transfer scores to species basis
	my $id_sp_ref = Pipeline::get_sp_from_faF($hmmFsR->{fasta});
	my %members; # all group member info
	my $n_multidom=0;
	my @ids = sort keys %hits; # queries are ids for seqs in fasta file
	my @allhrss;
	foreach my $id (@ids) {
		die "no species for id $id \n" unless (exists $id_sp_ref->{$id});
		my @ts = keys %{ $hits{$id}}; my $t = $ts[0]; # only 1 target id for self hmmscans
		$n_multidom++ if ($hits{$id}{$t}{rep} > 1); 
		my $ndom = $hits{$id}{$t}{rep};
		my $hrss = $hits{$id}{$t}{best_score} / $maxposssc;
		push @allhrss, $hrss;
		my $sp = $id_sp_ref->{$id};
		$members{$sp} = {id=>$id, hrss=>$hrss, ndom=>$ndom}; # also stored for each group -  a hash of species => id relscore ndom relscore_dev

	}
	my ($av, $min, $max, $sd) = _av_min_max_sd(@allhrss);
	foreach my $sp (keys %members) {
		if ($sd == 0) { # all scores identical (can occur in v rare cases)
			$members{$sp}{hrss_dev} = 0;
		} else {
			$members{$sp}{hrss_dev} = ($members{$sp}{hrss} - $av) / $sd;
		}

	}
	#-
	#+ --- put everything in output structure
	my %grpdetails;
	$grpdetails{maxposs_score} = $maxposssc;
	$grpdetails{mingrp_relscore} = $min;
	$grpdetails{avgrp_relscore} = $av;
	$grpdetails{maxgrp_relscore} = $max;
	$grpdetails{n_mem} = scalar keys %members;
	$grpdetails{n_multidom} = $n_multidom;
	$grpdetails{members} = \%members;
	#-
	# tidy up and return
	_remove_hmmfiles($hmmFsR, \@KEEPDEF); # files created by HMMER no longer needed, except specified
	return \%grpdetails;
}

sub get_alt_scores_groupdb{ 
# run a query fasta file versus group db and return hmrs scores for each query id
# (used to run hmmscan for alt orthos versus HMM0 grp dbs)
# my $altid_scoreR = HMMs::get_alt_scores_groupdb($grpid, $altfaF);
	my ($grpid, $altfaF) = @_;	
	my %altid_score; # holds all results
	if (-e $altfaF) {
		my $hmmFsR = Pipeline::group_fnames($grpid);
		my $hmmF = $hmmFsR->{hmm};
		my $althmmscanF = $hmmFsR->{althmmscan};
		my $althmmscandomF = $hmmFsR->{althmmscandom};
		_run_hmmscan($althmmscanF, $althmmscandomF, $hmmF, $altfaF);
		my %althmmFs = (hmmscan=>$althmmscanF, hmmscandom=>$althmmscandomF);
		my %althits = %{ get_hmm_scores( \%althmmFs) };
		# %althits is detailed info . Only need to return best score for each alt id %altid_score
		foreach my $id (keys %althits) {
			my @ts = keys %{ $althits{$id}}; my $t = $ts[0]; # only 1 target id for grp hmmscans
			$altid_score{$id} = $althits{$id}{$t}{best_score};
		}
		#-
		# tidy up and return
		_remove_hmmfiles($hmmFsR, \@KEEPDEF); # files created by HMMER no longer needed, except specified
	}
	return \%altid_score;
}

sub try_replacing_seq_in_group{
# usage: my ($isbetter, $resulttxt, $grpdetailsR) = HMMs::try_replacing_seq_in group($grpid, $oldgrpdetailsR, $repsp, $altfaF, $altidtotry);
# OR:    my ($isbetter, $resulttxt, $grpdetailsR) = HMMs::try_replacing_seq_in group($grpid, $oldgrpdetailsR, $repsp, $gbfaF, $genblastid);
	my ($grpid, $oldgrpdetailsR, $repsp, $repfaF, $id) = @_;
#+ get seq to try replacing with
	my ($idsR, $id_seqR) = Pipeline::read_fasta( $repfaF );
	die "id to try replacing  with for $grpid $repsp : $id - not found in $repfaF\n" unless (exists $id_seqR->{$id});
	my $repseqobj = $id_seqR->{$id};
	my $repdesc = $repseqobj->desc();# add sp to descriptor
	$repdesc ="species:$repsp " . $repdesc;
	$repseqobj->desc($repdesc); 
#-
#+ get group seqs
	my $hmmFs_ref = Pipeline::group_fnames($grpid);
	my $faF = $hmmFs_ref->{fasta};

	my $in = new Bio::SeqIO( -format => 'fasta', -file   => $faF );
	my %sp_seqobj;
	while (my $seqobj = $in->next_seq() ) {
		my $desc = $seqobj->desc();
		$desc =~ /species:([^ ]+) / || die "species info not found in $desc\n";
		die "multiple seq for sp $1 in group fasta $faF\n" if (exists $sp_seqobj{$1});
		$sp_seqobj{$1} = $seqobj;
	}
#-
#+ replace grp seq for sp
	change_hmm_files($grpid, "copy"); # keep copy of hmm files to revert if no improvement
	$sp_seqobj{$repsp} = $repseqobj;
	my $write_hitseq = new Bio::SeqIO( -format => 'fasta', -file   => ">$faF" );
	foreach my $sp (sort keys %sp_seqobj) {
		$write_hitseq->write_seq($sp_seqobj{$sp}); 
	}
#-
#+ redo HMMER and compare old and new scores
	my $grpdetailsR = HMMs::make_groupdb($grpid);
	my $newsphrss = $grpdetailsR->{members}{$repsp}{hrss};
	my $newscore = $newsphrss * $grpdetailsR->{maxposs_score};
	my $newavgrprelscore = $grpdetailsR->{avgrp_relscore};
	my $oldsphrss = 0; my $oldscore = 0; my $oldavgrprelscore = 0;
	if (defined $oldgrpdetailsR->{members}{$repsp}{hrss}) {
		$oldsphrss = $oldgrpdetailsR->{members}{$repsp}{hrss};
		$oldscore = $oldsphrss * $oldgrpdetailsR->{maxposs_score};
		$oldavgrprelscore = $oldgrpdetailsR->{avgrp_relscore};
	}
	my $scorecmptxt = sprintf "%s %.1f %s %.1f", "$repsp scores: old ", $oldscore, "new ", $newscore;
#-
#+ criteria for better: if hrss improved by 0.5 or if abs score improved by mindif without damaging grpav_relscore 
	my $mindif = $MINHRSSDIF * $oldgrpdetailsR->{maxposs_score};
	my $isbetter = ( (($newsphrss - $oldsphrss) > 0.5) || 
	                 ( (($newscore - $oldscore) > $mindif) && (($newavgrprelscore - $oldavgrprelscore) > -0.01) ) ); # rarely, changing pep can damage overall profile so this avoids this
#-
#+ depending on results, replace or revert
	if ($isbetter) {
		change_hmm_files($grpid, "delete_copy"); # keep new version only
	} else {
		$grpdetailsR = $oldgrpdetailsR; #reset to existing
		change_hmm_files($grpid, "revert"); # go back to previous version
	}
	return ($isbetter, $scorecmptxt, $grpdetailsR);
}


sub change_hmm_files {
# called	change_hmm_files($grpid, "copy"); or change_hmm_files($grpid, "revert");
	my ($grpid, $command) = @_;
	my $hmmFs_ref = Pipeline::group_fnames($grpid);
	my $hmmFscpy_ref = Pipeline::group_fnames("$grpid\.copy");	
	foreach my $k (keys %$hmmFs_ref) {
		my $current = $hmmFs_ref->{$k};
		next unless (-e $current); 
		my $copy = $hmmFscpy_ref->{$k};
		if ($command eq "copy") {
			`cp $current $copy`;
		} elsif ($command eq "revert") {
			`rm $current`;
			`cp $copy $current`;
			`rm $copy`;		
		} elsif ($command eq "delete_copy") {
			`rm $copy`;		
		} else {
			die "unrecognised command $command in change_hmm_files\n";
		}
	}
}

		
sub _remove_hmmfiles {
#
	my ($hmmFsR, $keeptheseR) = @_;
	my %keep;
	@keep{@$keeptheseR} = 1; # keep specified files per group
	foreach my $k (keys %$hmmFsR) {
		unless (exists $keep{$k}) {
			my $f = $hmmFsR->{$k};
			`rm -f $f`;
		}
	}
	unless (exists $keep{"hmm"}) {
		my $hmmF = $hmmFsR->{hmm};
		`rm -f $hmmF.*`; # remove binary files created by hmmpress
	}

}



sub make_db {
# make a HMMER db of set of groups with current hmm
# called my $dbmsg = HMMs::make_db($db, $dbgrpidsR); 
	my ($db, $dbgrpidsR) = @_; #  root name for all db files and set of grps in db
	#+ check whether db there
	my $fsthere = 1;
	my @dbFs;
	my @DBFEXTS = ("hmm", "hmm.h3f", "hmm.h3i", "hmm.h3m", "hmm.h3p"); # db file exts
	foreach my $ext (@DBFEXTS) {
		my $F = $HMMDIR . $db . "." . $ext;
		if (-e $F) {
			push @dbFs, $F;
		} else { 
			$fsthere = 0;
			last;
		}
	}
	#-
	my $dbmsg;
	if ($fsthere) {
		warn "Files for $db already present. Delete these if you want to recreate:\n";
		foreach my $dbF (@dbFs) {warn "$dbF\n";}
		$dbmsg = "HMMER database $db present";
	} else {
		chdir $HMMDIR;
		my $dbstoF = $db . ".sto";
		`rm $dbstoF` if (-e $dbstoF); 
		foreach my $grpid (sort keys %$dbgrpidsR) {
			my $hmmFsR = Pipeline::group_fnames($grpid);
			my $stoF = $hmmFsR->{sto};
			die "make_db: stockholm msa file $stoF not found\n" unless (-e $stoF);
			`cat $stoF >> $dbstoF`; # append stockholm format MSA
		}
		die "make_db: stockholm format msa file $dbstoF not found\n" unless (-e $dbstoF);
		my $dbhmmF = $dbstoF; $dbhmmF =~ s/\.sto/\.hmm/;
		`hmmbuild --amino --fragthresh 0 $dbhmmF $dbstoF`;
		die "make_db: No hmmbuild output file $dbhmmF created\n" unless (-e $dbhmmF);
		`hmmpress -f $dbhmmF`;
		my $nprof = keys %$dbgrpidsR;
		$dbmsg = "HMMER database $db created with $nprof profiles";
		chdir $TOPDIR;
	}
	return $dbmsg;
}

sub dbname_grpset {
# generate database name and set of grps to use according to supplied options
# called my ($db, $dbgrpidsR) = HMMs::dbname_grpset($grpid_statusR, $chmm, $ssopt);
	my ($grpid_statusR, $chmm, $ssopt) = @_;
	my $db;
	if ($chmm eq "HMM1") {
		$db = "HMM1";
		$ssopt = "all";
	} elsif ($chmm eq "final_hmm") {
		$db = "final_db";
		$db .= "_$ssopt" unless ($ssopt eq "all");		
	} else {
		die "dbname_grpset: unrecognised current hmm option $chmm\n";
	}
	my @checkgrpids;
	if ($ssopt eq "all") {
		@checkgrpids = keys %$grpid_statusR;
	} else { # $ssopt is interpreted as subset filename
		my ($ssgrpidsR, $ssdb) = Pipeline::read_subset($ssopt, $grpid_statusR); 
		@checkgrpids = @$ssgrpidsR;
	}
	#+ check that all grps are tested and find passsed 
	my @notdonegrpids;
	my %dbgrpids;
	foreach my $grpid (@checkgrpids) {
		if ($grpid_statusR->{$grpid}{"current_hmm"} ne $chmm) {
			push @notdonegrpids, $grpid;
		} else {
			if ($chmm eq "HMM1") {
				$dbgrpids{$grpid} = 1;
			} else { # other checks for inclusion in final HMM db
				my ($pass, $stattxt) = final_hmmstatus($MINRELSC, $grpid_statusR->{$grpid}{"current_hmm_details"});
				$dbgrpids{$grpid} = 1 if ($pass);
			}
		}
	}
	die "dbname_grpset: not all grpids for $db have expected hmm of $chmm . Following are incorrect stage :\n@notdonegrpids\n" if (@notdonegrpids > 0);
	#-
	return ($db, \%dbgrpids);
}
	

sub get_hmm_scores {
# called		my %hits = %{ get_hmm_scores($hmmFsR) };
# or            my $hits_ref = HMMs::get_hmm_scores($hmmFsR);
	my $hmmFsR = shift();
	my $wholetableR = parse_hmmscan($hmmFsR, "hmmscan");
	my %hits;
	my %multidom; # keep hits with with multiple domains, need to read in per_domain output
	foreach my $col (@$wholetableR) {
		my $q = $col->{query};
		my $t = $col->{target};	
		die "? >1 row for same query $q - target $t in " . $hmmFsR->{"hmmscan"} . " \n" if (exists $hits{$q}{$t});
		foreach my $k ("score_full", "score_1dom", "rep") {
			$hits{$q}{$t}{$k} = $col->{$k};
		}
		$multidom{$q}{$t} = $hits{$q}{$t}{"rep"} if ($hits{$q}{$t}{"rep"} > 1);
	}
	if (keys %multidom > 0) {	
		$wholetableR = HMMs::parse_hmmscan($hmmFsR, "hmmscandom");
		my %q_t_doms;
		foreach my $col (@$wholetableR) {
			my $q = $col->{query};
			my $t = $col->{target};	
			next unless (exists $multidom{$q}{$t});
			my $d = $col->{domnum};	# domain number
			my %dinfo;
			foreach my $k ("score_dom","bias_dom","from_hmm_coord","to_hmm_coord","from_ali_coord","to_ali_coord") {
				$dinfo{$k} = $col->{$k};
			}
			push @{ $q_t_doms{$q}{$t} }, \%dinfo;
		}
		foreach my $q (keys %multidom) {
			foreach my $t (keys %{ $multidom{$q}} ) {
				my $nf = scalar @{ $q_t_doms{$q}{$t} };
				die "found $nf doms for $q $t in " . $hmmFsR->{"hmmscandom"} . " expecting $multidom{$q}{$t} \n" unless 
				    ($nf == $multidom{$q}{$t});
				my ($mds, $ncompdoms) = _multidom_score( $q_t_doms{$q}{$t} );
				if ($ncompdoms == 1) { # only 1 dom found so just use single dom score
					$hits{$q}{$t}{"best_score"} = $hits{$q}{$t}{"score_1dom"};
				} elsif ($ncompdoms == $multidom{$q}{$t}) { # all doms found compatible so use full score
					$hits{$q}{$t}{"best_score"} = $hits{$q}{$t}{"score_full"};
				} else {
					$hits{$q}{$t}{"best_score"} = $mds;
				}
			}
		}	
	}
	foreach my $q (keys %hits) {
		foreach my $t (keys %{ $hits{$q}} ) {
			if (exists $hits{$q}{$t}{"best_score"}) {
# #+ for testing
#				warn "for $q $t full_score :". $hits{$q}{$t}{"score_full"} .
#				" multidom score :". $hits{$q}{$t}{"best_score"} ." \n";
# #-				
			} else {
				$hits{$q}{$t}{"best_score"} = $hits{$q}{$t}{"score_full"};
# #+ for testing
#				warn "for $q $t full_score :". $hits{$q}{$t}{"score_full"} ." \n";
# #-				
			}
		}
	}
	return \%hits;
}

sub _multidom_score {
	my $TOL=16; # tolerance - doms with overlaps > $TOL are eliminated as not compatible
# called 	my ($mds, $ncompdoms) = multidom_score( $q_t_doms{$q}{$t} );
	my $domsref = shift();
	my @doms = sort {$b->{score_dom} <=> $a->{score_dom}} @$domsref; # sort by descending scores
	# ----- find consistent multiple doms
	my @c_doms; # consistent doms
	push (@c_doms, $doms[0]);
	my %keep_dom;
	for (my $i=1; $i<@doms; $i++) { #  initialise with indeces of all other doms
		$keep_dom{$i}=1;
	}
	while (scalar(keys %keep_dom)>0) {
		foreach my $i (sort keys %keep_dom) {
			my $d1 = $doms[$i]; # check this one for _overlaps with all in c_doms
			foreach my $d2 (@c_doms) {
				if (Pipeline::overlap($d1->{from_hmm_coord}, $d1->{to_hmm_coord}, $d2->{from_hmm_coord}, $d2->{to_hmm_coord})>$TOL || 
				Pipeline::overlap( $d1->{from_ali_coord}, $d1->{to_ali_coord}, $d2->{from_ali_coord}, $d2->{to_ali_coord})>$TOL || 
				# check q and t order of doms are same 
				(($d2->{from_hmm_coord}-$d1->{from_hmm_coord})*($d2->{from_ali_coord}-$d1->{from_ali_coord})) < 0) {
					delete $keep_dom{$i};
					last;
				}
			}
			if (exists $keep_dom{$i}) {
				push (@c_doms, $doms[$i]); 
				delete $keep_dom{$i};
			}
		}	
	}
	my $mds = 0;
	my $ncompdoms = scalar @c_doms;
	foreach my $c_dom (@c_doms) {
		$mds += $c_dom->{score_dom} - $c_dom->{bias_dom};
	}
	return ($mds, $ncompdoms);
}

sub _av_min_max_sd {
	my @vals = @_;
	my $sum=0;
	my $min=$vals[0];
	my $max=$vals[0];
	foreach my $v (@vals) {
		$sum += $v;
		$min = $v if ($v < $min);
		$max= $v if ($v > $max);
	}
	my $n = scalar @vals;
	my $av = $sum / $n;
	my $ss=0;
	foreach my $v (@vals) {
		my $r = $v - $av;
		$ss += $r*$r;
	}
	my $sd = sqrt( $ss / ($n-1));
	return ($av, $min, $max, $sd);
}

sub parse_hmmscan {
# called: my ($coltitles_ref, $wholetableR) = HMMs::parse_hmmscan($hmmFsR, $option);
# $option indictaes per_hit or _per_domain hmmscan output
	my ($hmmFsR, $option) = @_;
	die "Unrecognised option $option\n" unless (exists $COLTITLES{$option});
	my $inF = $hmmFsR->{$option};
	die "File $inF does not exist in parse hmmscan\n" unless (-e $inF);
	open (IN, $inF) || die "Cannot open $inF for reading\n";
	my @coltitles = @{ $COLTITLES{$option}};
	my @wholetable;
	while (my $line = <IN>) {
		chomp $line;
		next if ($line =~ /^#/);
		my @a = split( /\s+/, $line);
		my %h;
		@h{ @coltitles } = @a;
		push @wholetable, \%h;
	}
	close IN;	
	return (\@coltitles, \@wholetable);
}

sub final_hmmstatus {
# called	my ($pass, $stattxt) = HMMs::final_hmmstatus($minrelsc, $grpdetailsR);
	my ($minrelsc, $grpdetailsR) = @_;
	my $pass=1;
	my $stattxt;
	my @missfinal;
	my @lowscore;
	my %members = %{ $grpdetailsR->{members}};
	foreach my $sp (@sps) {
		if (! exists $members{$sp}) {
			push @missfinal, $sp;
			$pass=0;
		} elsif ($members{$sp}{hrss} < $minrelsc) {
			my $lowtxt = sprintf "%s %.2f  ", $sp, $members{$sp}{hrss};
			push @lowscore, $lowtxt;
			$pass=0;
		}
	}
	if ($pass) {
		$stattxt = sprintf "%s %.2f", "pass:min_cutoff", $minrelsc;
	} else {
		my $missingtxt = ""; my $lowscoretxt = "";
		$missingtxt = "missing @missfinal" if (@missfinal>0);
		$lowscoretxt = "too low score @lowscore" if (@lowscore>0);
		$stattxt = "fail:[no missing allowed. min_score $minrelsc\]:" . join(":", ($missingtxt, $lowscoretxt));
	}
	return ($pass, $stattxt);
}

sub run_query_hmmscan {
# called my $hmmFsR = HMMs::run_query_hmmscan($qsp, $qF, $db);
# output put on separate dir (unlike hmmscan output of group members)
	my ($qsp, $qF, $dbname) = @_;
	my $db = $HMMDIR . $dbname . ".hmm";
	die "cannot find expected HMMER db file $db\n" unless (-e $db);
	my $hmmscanF = $QSCANDIR . "$qsp\_$dbname\.hmmscan.txt";
	my $hmmscandomF = $QSCANDIR . "$qsp\_$dbname\.hmmscan.dom.txt";
	if (-e $hmmscanF) {
		print "\n\n>>>>>>>>>>>>> hmm scan output file already exists :\n$hmmscanF\nUsing existing file\n";
		print "remove this file then re-run to create new version.<<<<<<<<<<<<<<<<<<<<<<<<\n\n";
	} else {
		_run_hmmscan($hmmscanF, $hmmscandomF, $db, $qF);
	}
	my %hmmFs = (hmmscan=>$hmmscanF, hmmscandom=>$hmmscandomF);
	return \%hmmFs;
}

sub _run_hmmscan {
	my ($hmmscanF, $hmmscandomF, $hmmF, $qF) = @_;
	`hmmscan  -E 1.e-7  --tblout $hmmscanF  --domtblout $hmmscandomF $hmmF $qF`;
	die "hmmscan failed. No output file $hmmscanF\n" unless (-e $hmmscanF);
	die "hmmscan failed. No output file $hmmscandomF\n" unless (-e $hmmscandomF);
}

sub get_hmmscan_hits_filename {
# called 	my $hmmscanthF = HMMs::get_hmmscan_hits_filename($qsp, $db, $outopt);
# or 	my $hmmscanthF = HMMs::get_hmmscan_hits_filename($qsp, $db);
	my $qsp = shift();
	my $db = shift();
	my $outopt = shift() || "topgrp"; #default
	my $hmmscanhitF = $QSCANDIR . "$qsp\_$db\.hmmscan.txt";
	my $substr = ".hmmscan_" . $outopt . "hits.";
	$hmmscanhitF =~ s/\.hmmscan\./$substr/;
	return $hmmscanhitF;
}

sub output_hmmscan_hits { 
# For a supplied query sp, output hmmscan hits organised by query then grp (target)
# used different modes - outputting only top hit for each query (used for optimisation of HMM1 with grass queries)
#                        outputting only top hit for each grp (used for finding top non-grass hit for each grp finalHMM and for non-member grass hits (ass peps)
#                        outputting all hits (used for finding other member grass hits of finalHMM to define supergrps)
# Currently no score cut-off here. This is implemented after all spp compiled together 
# called 	my $outputmsg = HMMs::output_hmmscan_hits($hitsR, $grpid_statusR, $qsp, $db, $outopt);
	my %OUTOPT_DESC = ("all"=>"all hits", "topq"=>"top hit for each query", "topgrp"=>"top hit for each group profile");
	my ($hitsR, $grpid_statusR, $qsp, $db, $outopt) = @_;
	die "output_hmmscan_hits: unrecognised output option $outopt\n" unless (exists $OUTOPT_DESC{$outopt});
	my %hits = %$hitsR;
	my %t_tophit;
	my %q_tophit;
	if ($outopt eq "topgrp") {
		foreach my $q (keys %hits) {
			foreach my $t (keys %{ $hits{$q}} ) {
				next if (exists $t_tophit{$t} && ($hits{$t_tophit{$t}}{$t}{"best_score"} > $hits{$q}{$t}{"best_score"}));
				$t_tophit{$t} = $q; 
			}
		}
	} elsif ($outopt eq "topq") {
		foreach my $q (keys %hits) {
			foreach my $t (keys %{ $hits{$q}} ) {
				next if (exists $q_tophit{$q} && ($hits{$q}{$q_tophit{$q}}{"best_score"} > $hits{$q}{$t}{"best_score"}));
				$q_tophit{$q} = $t; 
			}
		}
	}
	my $outF = get_hmmscan_hits_filename($qsp, $db, $outopt);
	open (OUT, ">$outF") || die "cannot open $outF\n";
# my @HMMSCANHITSCOLS = ("query", "grpid", "hit_score", "maxposs_score", "hit_relscore", "lowest_mem_relscore", "diff_lowest_hit_relscore");
	print OUT join("\t", @HMMSCANHITSCOLS) . "\n";
	my $nout=0;
	my %grpfound;
	foreach my $q (sort keys %hits) {
		foreach my $t (sort keys %{ $hits{$q}} ) {
			next if ($outopt eq "topgrp" && $t_tophit{$t} ne $q);
			next if ($outopt eq "topq" && $q_tophit{$q} ne $t);
			my $grpid = $t; $grpid =~ s/.msa//;
			unless (exists $grpid_statusR->{$grpid}) { # no output unless grpid exists
				warn "? profile id (target) $t not found ( $grpid )\n";
				next; 
			}
			my @row = ($q, $grpid, $hits{$q}{$t}{"best_score"});
			my $mps = $grpid_statusR->{$grpid}{"current_hmm_details"}{"maxposs_score"};
			push @row, $mps;
			my $relsc = $hits{$q}{$t}{"best_score"} / $mps;
			push @row, $relsc;
			my $lss = $grpid_statusR->{$grpid}{"current_hmm_details"}{"mingrp_relscore"};
			push @row, $lss;
			push @row, $lss - $relsc;
			print OUT join("\t", @row) . "\n";
			$grpfound{$grpid}++;
			$nout++;
		}
	}
	close OUT;
	my $nprof = keys %grpfound;
	return "Wrote $OUTOPT_DESC{$outopt} for $nout hits of $nprof profiles to $outF\n";
}


sub get_nonmem_hitscores {
	# callled my ($db, $dbgrpidsR, $sp_grpid_asspepR, $sp_grpid_bettermatchesR, $grpid_qgrpid_spR) = HMMs::get_nonmem_hitscores($grpid_statusR, $chmm, $hitsopt, $cutoff);
	my ($grpid_statusR, $chmm, $hitsopt, $cutoff) = @_;
	my ($db, $dbgrpidsR) = dbname_grpset($grpid_statusR, $chmm, "all");
	#+ following is to get relevant member details of grpids used to make db and put in convenient structures
	my %membid_grpid;
	my %sp_grpid_hrss; # to check for better scores
	my $isfinal = ($chmm eq "final_hmm"); # final hmm; all spp should have members
	foreach my $grpid (sort keys %$dbgrpidsR) {
		foreach my $sp (@sps) {
			my $member;
			if (exists $grpid_statusR->{$grpid}{"current_hmm_details"}{"members"}{$sp}) {
				$member = $grpid_statusR->{$grpid}{"current_hmm_details"}{"members"}{$sp};
				my $membid = $member->{"id"};	
				$sp_grpid_hrss{$sp}{$grpid} = $member->{"hrss"};
				$membid_grpid{$membid} = $grpid;
			} else {
				die "get_nonmem_hitscores: no $sp member for $grpid\n" if ($isfinal);
			}
		}
	}
	#+ read in all files and store ass and member matches separately
	my %sp_grpid_bettermatches; # all match qs that are not grp members and score better : diff rel scores 
	my %sp_grpid_asspep; # all match qs that are not grp members : rel scores
	my %grpid_qgrpid_sp; # all different grpids identified from match q member of different grp - for defining supergrps

	foreach my $qsp (@sps) {
		my $inF = HMMs::get_hmmscan_hits_filename($qsp, $db, $hitsopt);
		open (IN, $inF) || die "get_nonmem_hitscores: cannot open $inF\n";
		<IN>;
		while (<IN>) {
			chomp;
# my @HMMSCANHITSCOLS = ("query", "grpid", "hit_score", "maxposs_score", "hit_relscore", "lowestmem_relscore", "diff_lowestmem_hit_relscore");
			my ($q, $grpid, $ghscore, $maxpossscore, $ghrelscore, @a) = split("\t", $_);
			die "get_nonmem_hitscores: Unrecognised group id $grpid in $inF\n" unless (exists $grpid_statusR->{$grpid});
			next if ($ghrelscore < $cutoff);
			if (exists $membid_grpid{$q}) { # match is a grp member
				my $qgrpid = $membid_grpid{$q};
				$grpid_qgrpid_sp{$grpid}{$qgrpid}{$qsp} = $ghrelscore unless ($qgrpid eq $grpid); # store if match is member of diff grp
			} else { # match is not a grp member
				$sp_grpid_asspep{$qsp}{$grpid}{$q} = $ghrelscore; 
				my $diffscore;
				if (exists $sp_grpid_hrss{$qsp}{$grpid}) {
					$diffscore = $ghrelscore - $sp_grpid_hrss{$qsp}{$grpid};
				} else {
					$diffscore = $ghrelscore;
				}
				$sp_grpid_bettermatches{$qsp}{$grpid}{$q} = $diffscore if ($diffscore>$MINHRSSDIF);
			}
		}
		close IN;
	}
	return ($db, $dbgrpidsR, \%sp_grpid_asspep, \%sp_grpid_bettermatches, \%grpid_qgrpid_sp);
		
}

sub output_nonmem_scores {
# called $nrow = HMMs::output_nonmem_scores($ASSOUTF, \%sp_grpid_asspep, \@dbgrpids);
	my ($outF, $tableR, $grpidsR) = @_;
	open (OUT, ">$outF") || die "cannot open $outF\n";
	my %table = %$tableR; 
	my @h = ("non member peps", @sps); # header row
	print OUT join("\t",@h) . "\n";
	my $n=0;
	foreach my $grpid (sort keys %$grpidsR) {
		my @r = ($grpid); #  start row
		my $hitinrow=0;
		foreach my $sp (@sps) {
			if (exists $table{$sp}{$grpid}) {
				my %q_s = %{ $table{$sp}{$grpid} };
				my @qs = sort {$q_s{$b} <=> $q_s{$a}} keys %q_s; # sort by descending score
				my @qstxts;
				foreach my $q (@qs) {
					my $qstxt = sprintf "%s %.3f",  $q, $q_s{$q};
					push @qstxts, $qstxt;
				}
				push @r, join(";", @qstxts);
				$hitinrow=1;
			} else {
				push @r, "";
			}
		}
		if ($hitinrow) {
			print OUT join("\t",@r) . "\n";
			$n++;
		}
	}
	close OUT;
	return $n;
}

sub input_nonmem_scores {
# called 	$grpid_sp_scoresR = HMMs::input_nonmem_scores($inF);
	my %grpid_sp_scores;
	my $inF = shift();
	open (IN, $inF);
	my $h = <IN>;
	while (<IN>) {
		chomp;
		my ($grpid, @persp) = split("\t", $_);
		my $i=0;
		foreach my $sp (@sps) {
			if (defined $persp[$i] && length $persp[$i] > 0) {
				my @qstxts = split (";", $persp[$i]);
				my %q_s;
				foreach my $qstxt (@qstxts) {
					my ($k, $v) = split (" ", $qstxt);
					$q_s{$k} = $v;
				}
				$grpid_sp_scores{$grpid}{$sp} = \%q_s;
			}
			$i++;
		}
	}
	close IN;
	return \%grpid_sp_scores;
}



1;
