##################################################
# universal_grass_peps pipeline v1.4             #
# Rowan Mitchell 2020-2024                       #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################

package Pipeline;
use strict;
use warnings;
require Exporter;

#+ ++++++++++++++++++++++++  exported
our @ISA = qw(Exporter);
our @EXPORT = qw($TOPDIR $BLASTPDIR $GMDIR $HMMDIR $GFFDIR $KEEPGBLASTDIR $QSCANDIR 
                 $REDSEQF @REFSPS %REFSP_SHTN $CLUSTFEXT $ZMOSF $ORTHOTABLEF $HMM1HITSF $GMSTATF $TODOTXT $COMPLETETXT
                 $SELFBLASTPLIM $NG_MAX_MISSGM $MINHRSSDIF 
                 @FSSTATFS $FSSTATFROOT $FSSTATFALL @FSSTAT_FIELDS @GRPDETSUM $OUTORTHOEXT @ORTHKEYS $MINRELSC $DEVLIM
                 @sps %sp_stat %refsp_clust_peps %refsp_pep_clust %refsp_clust_allorth %pepset %orth_topclust %grp_table);


#+ constants
#+ dirs
our $TOPDIR = "/home/data/tropical_forage_grasses/grass_core/final/";
our $BLASTPDIR = $TOPDIR . "blastp/";
our $GMDIR = $TOPDIR . "EnsemblPlants/pep_gene_models/";
our $FAHMM0DIR = $TOPDIR . "EnsemblPlants/orthologs/fasta_for_HMM0/";
our $FAALTDIR = $TOPDIR . "EnsemblPlants/fasta_for_altpeps/";
our $HMMDIR = $TOPDIR . "hmms/";
our $GFFDIR = $TOPDIR . "EnsemblPlants/gff/";
our $KEEPGBLASTDIR = $TOPDIR . "genblast/";
our $QSCANDIR = $HMMDIR . "query_hmmscans/";
our $NONGRASSDIR = $GMDIR . "nongrass/";
#- 
#+ filenames and filename parts
our $REDSEQF = $GMDIR . "redundant_ids.txt"; # record of seqs that have been removed (used for lookup in ortho comparison)
our $ZMOSF =  $BLASTPDIR . "Zea_mays-Os.hits.txt";
our $CLUSTFEXT = "-clusters.txt";
our $ORTHOTABLEF = $TOPDIR . "EnsemblPlants/orthologs/ortho_table.txt";
our $OUTORTHOEXT = "_allorths.txt";
our $HMM1HITSF = "HMM1_bettermatches.txt"; # for other peps that would apparently be better matches than existing member at HMM1 stage
our $NONGRASSINFO = $NONGRASSDIR . "nongrass_spp.txt";
#-
#+ for Ortho calc
our @REFSPS = ("Oryza_sativa", "Zea_mays");
our %REFSP_SHTN = ("Oryza_sativa"=>"Os", "Zea_mays"=>"Zm");
our @ORTHKEYS = ("q_pcid", "confidence", "t_pcid"); 
#-
#+ for status files
our $GMSTATF = $TOPDIR . "grass_sp_status.txt";
our $TODOTXT = "to_do";
our $COMPLETETXT = "complete";
#our @GMSTAT_FIELDS = ("taxon", "ori_pep_fasta_file", "genome_fasta_file", "gff3_file", "gmstep_completed");
our @GMSTAT_FIELDS = ("taxon", "ori_pep_fasta_file", "genome_files_start", "genome_files_chrs", "gff3_file", "gmstep_completed");	

our $FSSTATFROOT = "complete-group-steps.stat"; # Final group steps status files
our @FSSTATFS;
our $NSTATFS = 12; 
for (my $i=0; $i<$NSTATFS; $i++) {
	$FSSTATFS[$i] = sprintf "%s%s%02u%s", $TOPDIR, $FSSTATFROOT, $i+1, ".txt"; 
}
our $FSSTATFALL = $TOPDIR . $FSSTATFROOT . ".txt";
# status file col titles : "group_id", @FSSTAT_FIELDS, @GRPDETSUM, @sps;
our @FSTEPS = ("start", "HMM1", "HMM2", "genblast");
our %FSTEP_CHMM = ("start"=>"none", "HMM1"=>"HMM1", "HMM2"=>"HMM2", "genblast"=>"final_hmm");
our @FSSTAT_FIELDS = ("n_mem_ortho", "step_complete", "current_hmm");
our @GRPDETSUM = ("maxposs_score", "mingrp_relscore", "avgrp_relscore", "maxgrp_relscore", "n_mem", "n_multidom");  # group summary fields
#-
#+ parameters
our $SELFBLASTPLIM = 90.; # %id >= this for self hit to clustered
our $NG_MAX_MISSGM = 4; # max num missing grass spp for grp to be tried for refinement
our $MINHRSSDIF = 0.01; # min diff in HRSS profile score for it to be classed as better
our $MINRELSC = 0.65; # min rel score for final HMM group to pass as highly conserved
our $DEVLIM = -2.0; # number of SDs for member to be classed as an outlier and to try genblast
#-
# package variables exported
our %refsp_clust_peps; # clustering of peps into groups from self blastp clust -> [@peps]
our %refsp_pep_clust; # clustering of peps into groups from self blastp	pep -> clust
our %refsp_clust_allorth; # all orthologs for each refsp clust
our %pepset; # all seed peps fron ref spp
our %orth_topclust; # top matching clust for each orth
our %grp_table ; # to hold final output of orthos
our @sps;
our %sp_stat;
#- -------------------------------

# package variables not exported
our @header;

sub read_gmsteps_status {	# 
# callled: my $gmfinished = Pipeline::read_gmsteps_status();	
	die "file $GMSTATF does not exist\n" unless (-e $GMSTATF);
	open (IN, $GMSTATF)|| die "Cannot open input file $GMSTATF\n";
	my $line = <IN>;
	chomp $line;
	while ($line =~ /^#/) {
		push @header, $line;
		$line = <IN>;
		chomp $line;
	} # last $line is @GMSTAT_FIELDS
	my $gmfinished = 1; # flag for steps complete for all spp;
	@sps = ();
	while (<IN>) {
		chomp;
		my ($sp, @a) = split ("\t", $_);
		push @sps, $sp;
		my %h;
		@h{@GMSTAT_FIELDS} = @a;
		$gmfinished =($gmfinished && ($h{"gmstep_completed"} eq $COMPLETETXT));
		$sp_stat{ $sp } = \%h;
	}
	close IN;
	return $gmfinished;
}

sub gmfasta_filename {
# called 	my $faF = gmfasta_filename( $sp ); # gene model pep fasta filename
	my $sp = shift();
	read_gmsteps_status() unless (exists $sp_stat{$sp} );
	my $pF = $sp_stat{$sp}{"ori_pep_fasta_file"};
	$pF =~ s/\.all\./\.nr\./ || die "substituting \'nr\' for \'all\' failed for $pF";
	my $faF = $GMDIR . $pF;
	return $faF;
}


sub update_gmsteps_status {	# 
# called: 	update_gmsteps_status()
	open (OUT, ">$GMSTATF")|| die "Cannot open output file $GMSTATF\n";
	foreach my $hl (@header) {
		print OUT $hl . "\n";
	}
	print OUT join("\t", ("species", @GMSTAT_FIELDS)) . "\n";
	foreach my $sp (@sps) {
		my %h = %{ $sp_stat{ $sp }};
		print OUT join("\t", ($sp, @h{@GMSTAT_FIELDS})) . "\n";
	}
	close OUT;
}

sub read_fasta {
# called 	my ($idsR, $id_seqR) = Pipeline::read_fasta( $faF );
	my $faF = shift();
	die "File $faF not found\n" unless (-e $faF);
	my $in = new Bio::SeqIO( -format => 'fasta', -file   => $faF );
	my @ids;
	my %id_seqobj;
	while (my $seqobj = $in->next_seq() ) {
		my $id = $seqobj->display_id();
		push @ids, $id;
		$id_seqobj{$id} =  $seqobj;
	}
	return (\@ids, \%id_seqobj);
}

sub read_fasta_lens {
# called 	my ($idsR, $id_lenR) = Pipeline::read_fasta_lens( $faF );
	my ($idsR, $id_seqR) = read_fasta( shift() );
	my %id_len;
	my %id_seqobj = %{$id_seqR };
	foreach my $id (@$idsR) {
		$id_len {$id} = $id_seqobj{ $id }->length();
	}
	return ($idsR, \%id_len);
}

sub write_fasta {
# called 		my $nseq = Pipeline::write_fasta($faF, \%id_seqobj );
	my ($faF, $id_seqobjR ) = @_;
	my %id_seqobj = %{ $id_seqobjR};
	my $n = scalar keys %id_seqobj;
	if ($n>0) {
		my $write_hitseq = new Bio::SeqIO( -format => 'fasta', -file   => ">$faF" );
		foreach my $id (sort keys %id_seqobj) {
			$write_hitseq->write_seq( $id_seqobj{$id} ); 
		}
	}
	return $n;
}

sub append_fasta {
# called 		my $nseq = Pipeline::append_fasta($faF, \%id_seqobj );
	my ($faF, $id_seqobjR ) = @_;
	my %id_seqobj = %{ $id_seqobjR};
	my $n = scalar keys %id_seqobj;
	if ($n>0) {
		my $write_hitseq = new Bio::SeqIO( -format => 'fasta', -file   => ">>$faF" );
		foreach my $id (sort keys %id_seqobj) {
			$write_hitseq->write_seq( $id_seqobj{$id} ); 
		}
	}
	return $n;
}

sub add_seq_to_grpfasta {
# called		Pipeline::add_seq_to_grpfasta($sp, $faF, $id_seqobjR->{$q} );
	my ($sp, $faF, $seqobj ) = @_;
	my $write_hitseq = new Bio::SeqIO( -format => 'fasta', -file   => ">>$faF" );
	my $desc = $seqobj->desc();	# add sp info to descriptor
	$desc ="species:$sp " . $desc;
	$seqobj->desc($desc);
	$write_hitseq->write_seq($seqobj); 
	return;
}

sub group_fasta_fnames {
# called	my $faFsR = Pipeline::group_fasta_fnames($grpid); 
# [could make alt and hmm1 alt the same so they overwrite and save disk space]
	my $grpid = shift();
	my $hmm0faF = $FAHMM0DIR . $grpid . ".HMM0.fa";
	my $altfaF = $FAHMM0DIR . $grpid . "_alt.HMM0.fa";
	my $hmm1altfaF = $FAALTDIR . $grpid . "_alt.fa";
	my %faFs = (HMM0=>$hmm0faF, HMM0alt=>$altfaF, HMM1alt=>$hmm1altfaF);
	return \%faFs;
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#vvvvvvvvvv  final steps status file routines
# status file cols : "group_id", "n mem (ortho_table)", @FINALSTEPS, "current_hmm", @GRPDETSUM, @sps;

sub create_fsteps_status_file {
#called Pipeline::create_fsteps_status_file();
# after completion of all gene model steps, creates empty status file 
# not populated here @grpdetsum, @grpdetsps;
	open (OUT, ">$FSSTATFALL");
	print OUT join("\t", ("group_id", @FSSTAT_FIELDS, @GRPDETSUM, @sps)) . "\n";
	my $step = $FSTEPS[0];
	my $chmm = $FSTEP_CHMM{$step};
	foreach my $grpid (sort keys %grp_table) {
		my $nsp = scalar keys %{ $grp_table{$grpid}};
		print OUT join("\t", ($grpid, $nsp, $step, $chmm)) . "\n";
	}
	close OUT;
	return;
}

		
sub update_fsteps_status {
# called		Pipeline::update_fsteps_status(\%grpid_status, $fsstatF);
	my ($grpid_statusR, $fsstatF) = @_;
	open (OUT, ">$fsstatF") || die "cannot open $fsstatF\n";
	print OUT join("\t", ("group_id", @FSSTAT_FIELDS, @GRPDETSUM, @sps)) . "\n";
	my %grpid_status = %$grpid_statusR;
	foreach my $grpid (sort keys %grpid_status) {
		my @row = ($grpid);
		my $step = $grpid_status{ $grpid}{"step_complete"};
		die "step_complete for grpid $grpid not recognised: $step\n" unless (exists $FSTEP_CHMM{$step});
		$grpid_status{ $grpid}{"current_hmm"} = $FSTEP_CHMM{$step}; # only field that is set here
		foreach my $f (@FSSTAT_FIELDS) {
			push @row, $grpid_status{ $grpid}{$f};
		}
		unless ($step eq $FSTEPS[0]) { # has hmm fields
			if (exists $grpid_status{ $grpid}{"current_hmm_details"}) { 
				my @gds = @{ _group_details_texts($grpid_status{ $grpid}{"current_hmm_details"})};
				@row = (@row, @gds); # append row with grp hmm details
			} else {
				my $chmm = $grpid_status{ $grpid}{"current_hmm"};
				warn "expected details of $chmm hmm not supplied for $grpid\n";
			}				
		}
		print OUT join("\t", @row). "\n";
	}
	close OUT;
	return;
}

sub get_fsteps_status {
# called my ($grpid_statusR, $fsstatF) = Pipeline::get_fsteps_status($mode);
	my $mode = shift();
	my $fsstatF;
	if ($mode eq "all") {
		$fsstatF = $FSSTATFALL;
	} else {
		die "Unrecognised set $mode - should be > 0 <= $NSTATFS\n" unless
		  ($mode>0 && $mode<=$NSTATFS);
		$fsstatF = $FSSTATFS[ $mode-1 ];
	}	
	die "file $fsstatF not there\n" unless (-e $fsstatF);
	open (IN, $fsstatF) || die "cannot open $fsstatF\n";
	my %grpid_status;
	<IN>; # col title row
	while (<IN>) {
		chomp;
		my ($grpid, @a) = split("\t", $_); # 
		my %stat;
		for (my $i=0; $i<@FSSTAT_FIELDS; $i++) {
			$stat{ $FSSTAT_FIELDS[$i] } = $a[$i];
		}
		unless ($stat{"step_complete"} eq $FSTEPS[0]) { # has hmm fields		
			my $i = @FSSTAT_FIELDS; 
			my $j = scalar @a - 1;
			my @gds = @a[$i..$j];
			$stat{"current_hmm_details"} = _parse_group_details(\@gds);
		}
		$grpid_status{ $grpid} = \%stat;
	}
	close IN;
	return (\%grpid_status, $fsstatF);
}

#+ routines to put and retrieve group details in texts for status file
sub _group_details_texts {
	my $grpdetailsR = shift();
	die "species list not initialised _group_details_texts\n" unless (@sps>0);
	my %gd = %{ $grpdetailsR};
	my @a = @gd{@GRPDETSUM};
	my %m = %{ $gd{members} }; # members
	for (my $i=0; $i<@sps; $i++) {
		my $j = $i + @GRPDETSUM;
		if (exists $m{$sps[$i]}) {
			my %me = %{ $m{$sps[$i]} };
			warn "? undefined grp details for sp $sps[$i] \n" unless (defined $me{"ndom"});
			$a[$j] = sprintf "%s %.3f %.1f %.0f", @me{("id", "hrss", "hrss_dev", "ndom")};
		} else {
			$a[$j] = "";
		}
	}
	return \@a;
}

sub _parse_group_details {
	my @a = @{ shift()};
	die "species list not initialised _parse_group_details\n" unless (@sps>0);
	my %gd;
	for (my $i=0; $i<@GRPDETSUM; $i++) {
		$gd{$GRPDETSUM[$i]} = $a[$i];
	}
	my %m; # members
	for (my $i=0; $i<@sps; $i++) {
		my $j = $i + @GRPDETSUM;
		next unless (defined $a[$j] && $a[$j] =~ /\w/);
		my ($id, $hrss, $hrss_dev, $ndom) =  split(" ", $a[$j]);
		$m{$sps[$i]} = {id=>$id, hrss=>$hrss, hrss_dev=>$hrss_dev, ndom=>$ndom};				
	}
	$gd{members} = \%m;
	return \%gd;
}
#-

#^^^^^^^  final steps status file routines
#-----------------------------------------------------------------------------------

sub read_subset {
# called:	my ($ssgrpidsR, $ssdb) = Pipeline::read_subset($subsetF, $grpid_statusR); 
	my ($subsetF, $grpid_statusR) = @_;
	die "file $subsetF not there\n" unless (-e $subsetF);
	open (IN, $subsetF) || die "cannot open $subsetF\n";
	my %grpids;
	while (<IN>) {
		next if (/^#/);
		chomp;
		my ($grpid, @a) = split("\t", $_); # 
		$grpids{ $grpid} = 1;
	}
	close IN;
	my @ssgrpids;
	foreach my $ssgrpid (keys %grpids) {
		if (exists $grpid_statusR->{ $ssgrpid}) {
			push @ssgrpids, $ssgrpid;
		} else {
			warn "id $ssgrpid from $subsetF does not exist\n";
		}
	}
	my $ssdb = $subsetF; $ssdb =~ s/\.txt//; # root name for db
	return (\@ssgrpids, $ssdb);
}
	
sub group_fnames {
# called	my $grpFsR = Pipeline::group_fnames($grpid);	
# or	my $grpFsR = Pipeline::group_fnames($grpid, $ext);	
	my $grpid = shift();
	my $useext = (my $ext = shift());
	my $root;
	if ($useext) {
		$root = $grpid . ".$ext";
	} else {
		$root = $grpid;
	}
	my $faF = $HMMDIR . $root . ".fa";
	my $mlaF = $faF; $mlaF =~ s/\.fa/\.msa\.fa/;
	my $hmmF = $faF; $hmmF =~ s/\.fa/\.hmm/;
	my $stoF = $faF; $stoF =~ s/\.fa/\.sto/;
	my $hmmscanF = $faF; $hmmscanF =~ s/\.fa/\.hmmscan\.txt/;
	my $hmmscandomF = $faF; $hmmscandomF =~ s/\.fa/\.hmmscan\.dom\.txt/;
	my $confaF = $hmmF . "_consensus.fa";
	my $conhmmscanF = $confaF; $conhmmscanF =~ s/\.fa/\.hmmscan\.txt/;
	my $conhmmscandomF = $confaF; $conhmmscandomF =~ s/\.fa/\.hmmscandom\.txt/;
	my $altfaF = $HMMDIR . $root . "_alt.fa";
	my $althmmscanF = $altfaF; $althmmscanF =~ s/\.fa/\.hmmscan\.txt/;
	my $althmmscandomF = $altfaF; $althmmscandomF =~ s/\.fa/\.hmmscandom\.txt/;
	my %grpFs = (fasta=>$faF, mla=>$mlaF, hmm=>$hmmF, sto=>$stoF, hmmscan=>$hmmscanF, hmmscandom=>$hmmscandomF, 
	consensus=>$confaF, conhmmscan=>$conhmmscanF, conhmmscandom=>$conhmmscandomF,
	alt=>$altfaF, althmmscan=>$althmmscanF, althmmscandom=>$althmmscandomF);
	return \%grpFs;
}

sub get_sp_from_faF {
# called: 	my $hitid_spR = Pipeline::get_sp_from_faF($inF);
# assigns each id in a fasta file to a grass species from descriptor info
	my $faF = shift();
	die "File $faF not found\n" unless (-e $faF);
	my $in = new Bio::SeqIO( -format => 'fasta', -file   => $faF );
	my %id_sp;
	while (my $seqobj = $in->next_seq() ) {
		my $id = $seqobj->display_id();
		my $desc = $seqobj->desc();
		$desc =~ /species:([A-Za-z0-9_\-]+) / || die "species info not found in $desc\n";
		$id_sp{$id} = $1;
	}
	return \%id_sp;
}

sub write_fsteps_subset { # convenient to have file just of subset
# called 	my $resF = Pipeline::write_fsteps_subset($setgrpidsR, \%grpid_status, $subsetF); # convenient to have file just of subset
	my ($setgrpidsR, $grpid_statusR, $subsetF) = @_;
	my $resF = $TOPDIR . $subsetF; $resF =~ s/\.txt/\.results\.txt/;
	open (OUT, ">$resF") || die "cannot open $resF\n";
	print OUT join("\t", ("group_id", @FSSTAT_FIELDS)) . "\n";
	my %grpid_status = %$grpid_statusR;
	foreach my $grpid (sort @$setgrpidsR) {
		my @row = ($grpid);
		if (exists $grpid_status{ $grpid}) {
			foreach my $f (@FSSTAT_FIELDS) {
				push @row, $grpid_status{ $grpid}{$f};
			}
		} else {
			push @row, "NOT IN SET OF GROUPS THAT PASSED GM STEPS";
		}
		print OUT join("\t", @row). "\n";
	}
	close OUT;
	return $resF;
}

sub nongrass_info {
# called 	my ($qspsR, $ngsp_infoR) = Pipeline::nongrass_info();
	open (NGF, $NONGRASSINFO) || die "cannot open $NONGRASSINFO\n";
	my @sps;
	my %ngsp_info;
	while (<NGF>) {
		chomp;
		my ($sp, $pf, $tx) = split("\t", $_);
		push @sps, $sp;
		$ngsp_info{$sp}{pep_fasta} = $NONGRASSDIR . $pf . ".pep.all.fa";
		$ngsp_info{$sp}{taxon} = $tx;
	}
	return (\@sps, \%ngsp_info);
	close NGF;
}

sub overlap {
	my ($a1, $a2, $b1, $b2) = @_;
	my $ol;
	if ($b1 > $a1) { # b after a
		$ol = $a2 - $b1;
	} elsif ($b1 < $a1) { # a after b
		$ol = $b2 - $a1;
	} else { # a and b start at same pos
		if ($a2 >= $b2) {
			$ol = $b2 - $a1;
		} else {
			$ol = $a2 - $a1;
		}
	}
	$ol = 0 if ($ol < 0);
	return $ol;
}

1;
