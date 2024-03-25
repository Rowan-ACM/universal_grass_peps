##################################################
# universal_grass_peps pipeline v1.4             #
# Rowan Mitchell 2020-2024                       #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################


package Genblast;

use strict;
use warnings;
use Pipeline;



my $GENBLASTDIR = "/home/data/bioinf_resources/programming_tools/genBlastA/genBlast_v138_linux_x86_64/";
my $GENBLASTSUBDIR = "grassgenes/";
my %GBLASTOFS_EXTS = ("pro" => ".gblast.tophit.pro", "cdna" => "gblast.tophit.DNA", "gff" => "gblast.tophit.gff", "info" => "gblast.tophit");
my %GFFPRIMARYTAGS = ("ensembl"=>{"cds"=>"mRNA", "exon"=>"CDS"},
                      "genblast"=>{"cds"=>"transcript", "exon"=>"coding_exon"});



sub run_genblast {
# called	my ($failed, $reporttxt, $gbpepfaF, $gbpepid) = Genblast::run_genblast( $sp, $grpid);
	my ($sp_to_redo, $grpid) = @_;
	my $hmmFsR = Pipeline::group_fnames($grpid);
	my $faF = $hmmFsR->{fasta};
	my $confaF = $hmmFsR->{consensus};
	die "File $faF not found\n" unless (-e $faF);
	die "File $confaF not found\n" unless (-e $confaF);
	chdir $GENBLASTDIR;
	my $qF = $GENBLASTSUBDIR . $grpid . "_consensus.fa";
	my $infa = new Bio::SeqIO( -format => 'fasta', -file   => $confaF );
	my $seqobj= $infa->next_seq();
	my $qid = $seqobj->display_id();# id of > 30 chars causes error in genblast so use in desc only
	$seqobj->desc("original id: $qid");
	$seqobj->display_id("hmm_consensus");
	# write to query fa file
	my $outfa = new Bio::SeqIO( -format => 'fasta', -file   => ">$qF" );
	$outfa->write_seq($seqobj); 
	#+ outputs
	my $failed = 1; # a flag for no hits found by genblast
	my $reporttxt; # text of outcome
	my $gbpepid; # id of tophit pep in fasta file
	my $gbpepfaF; #  fasta file with hit copied from genblast
	#- 
	#+ tophit vars to find from loop
	my $tophitscore = 0.0; 
	my $tophitFsR = []; #  genblast files for tophit to rename and keep
	my $tophitseqobj;
	my $topchr;
	#-
	#+ for Ensembl gene model details
	my $gffF = $GFFDIR . $sp_stat{$sp_to_redo}{"gff3_file"};
	my %chr_existinggms;
	#-
	# multiple chrs loop -------------
	foreach my $chr (split (";", $sp_stat{$sp_to_redo}{"genome_files_chrs"})) {
		my $gFshort = $sp_stat{$sp_to_redo}{"genome_files_start"} . $chr . ".fa";
		my $gF = $GENBLASTSUBDIR . $gFshort;
		die "File $gF not found\n" unless (-e $gF);
		my $oF = $GENBLASTSUBDIR . $grpid . "_" . $gFshort; # output files stem
		`rm -f $oF*`; # remove any old genblast output files for this query-target combination
		`./genblast_v138_linux_x86_64  -p genblastg -q $qF -t $gF -o $oF -v 2 -h 0 -j 3 -r 1  -norepair -gff -cdna -pro -pid`;
		my @oFs  = `ls $oF*`;
		my @nFs;
		foreach my $oldname (@oFs) {  # tidy up names
			$oldname =~ s/\s//g;
			my $newname = $oldname;
			$newname =~ s/_1\.1c_2\.3_s2_tdshift0_tddis0_tcls3\.0_m2_score_i0_d16_0/\.gblast/; 
			`mv $oldname $newname`;
			push @nFs, $newname;
		} 
		my $gbiF = $oF . ".gblast";
		die "Die: no main genblast output file $gbiF found\n" unless (-e $gbiF);
		_correct_giant_file( $gbiF ); # correct for genblast bug
		my $proF = $oF . ".gblast.pro";
		my $score=0;
		if (-e $proF) { # no hit unless genblast protein output file exists
			my $hitseqobj;
			# read in genblast hit
			my $in = new Bio::SeqIO( -format => 'fasta', -file   => $proF );
			if ($hitseqobj= $in->next_seq()) {; #failed if genblast protein output file empty
				$score = _parse_genblast_score( $gbiF);
				if ($score > $tophitscore) {
					unless (exists $chr_existinggms{$chr}) { # checks if gene details already read in
						my $gffppF = $gffF . "\.$chr\_pp.txt";
						$chr_existinggms{$chr} = get_existing_gms( $gffppF ); # read in preprocessed file of all gene models for this sp and chr
					}
					my $gffF = $oF . ".gblast.gff";
					my $genedet = genedet_genblast($gffF );
					if (exists $chr_existinggms{$chr}{$genedet}) {
						$reporttxt .= "\t$gFshort $score BUT identical to existing gene model " . $chr_existinggms{$chr}{$genedet} . "\n";
					} else {
						$failed = 0; # successfully obtained hit with +ve score so failed set false
						$tophitscore = $score; 
						$tophitFsR = [ @nFs]; # top hit files 
						$tophitseqobj = $hitseqobj;
						$topchr = $chr;
					}
				}
			}
		}
		if ($score > 0) {
			$reporttxt .= "\t$gFshort $score\n";	
		} else {
			$reporttxt .= "\t$gFshort no hit\n";
		}
	}	# end multiple chrs loop ---------------
	#+ rename tophit files to potentially keep
	foreach my $thF (@$tophitFsR) {
		my $newname = $thF;
		$newname =~ s/\.gblast/\.gblast\.tophit/;
		`mv $thF $newname`;
	}
	#-
	chdir $TOPDIR;
	unless ($failed) { 
		# change features of genblast hit 
		my $seqtxt = $tophitseqobj->seq();# remove * from seq (generates warning)
		$seqtxt =~ s/\*//;
		$tophitseqobj->seq($seqtxt);
		$gbpepid = "genblast_" . $grpid . "_" . $sp_to_redo . "_" . $topchr;
		$tophitseqobj->display_id($gbpepid);
		# write to replacement fa file
		$gbpepfaF = $HMMDIR . $gbpepid . ".tmp.fa";
		my $write_hitseq = new Bio::SeqIO( -format => 'fasta', -file   => ">$gbpepfaF" );
		$write_hitseq->write_seq($tophitseqobj); 
	}
		
	return ($failed, $reporttxt, $gbpepfaF, $gbpepid); 
}

sub _parse_genblast_score {
# called 			my $score = _parse_genblast_score( $scoreF);
	my $scoreF = shift();
	die "Die: Genblast score file $scoreF not found\n" unless (-e $scoreF);
	open (IN, $scoreF) || die "cannot open $scoreF \n";
	my $l;
	my $score;
	while ($l = <IN>) {
		if ($l =~ /Gene:.*rank:1.score:([-\d\.e]+).PID/) {
			$score = $1;
			last;
		}
	}
	close IN;
	unless (defined $score) {
		warn "could not parse score from $scoreF\n";
		$score = 0;
	}
	return $score;
}

sub keep_genblast_files {
# called Genblast::keep_genblast_files($grpid); # keep successful genblast output
	my $grpid = shift();
	foreach my $k (keys %GBLASTOFS_EXTS) {
		my $gbFsearch = $GENBLASTDIR . $GENBLASTSUBDIR . $grpid . "*" . $GBLASTOFS_EXTS{$k};
		my @gbFs = `ls $gbFsearch`;
		my $nfound = @gbFs;
		unless ($nfound == 1) {
			warn "found $nfound files with $gbFsearch - not keeping this genblast $GBLASTOFS_EXTS{$k} file\n";
			next;
		}
		my $gbF = $gbFs[0];
		$gbF =~ s/\s//g;
		my $oldfile = $gbF;
		unless (-e $oldfile) {
			warn "$oldfile not found - not keeping this genblast $GBLASTOFS_EXTS{$k} file\n";
			next;
		}
		$gbF =~ s/$GENBLASTDIR//;
		$gbF =~ s/$GENBLASTSUBDIR//;
		my $newfile = $KEEPGBLASTDIR . $gbF;
		`cp $oldfile $newfile`;
	}
	return;
}

sub _correct_giant_file { # check for odd huge .gblast file produced by bug in genblast
	my $gbiF = shift();
	my $lsl = `ls -s $gbiF`;
	my ($size, $fname) = split (" ", $lsl);
	if ($size > 10000) { 		
		warn "correcting huge file $gbiF\n";		
		my $tmpF = $gbiF . ".tmp";
		`rm -f $tmpF` if (-e $tmpF);
		open (IN, $gbiF);
		open (OUT, ">$tmpF");
		while (<IN>) {
			my $l = length($_);
			print OUT $_ if ($l<10000);
		}
		close IN;
		close OUT;
		`rm -f $gbiF`;
		`mv $tmpF $gbiF`;
	}
}
	
sub remove_genblast_files {
# called 		Genblast::remove_genblast_files($grpid); # remove all output on genblast dir
	my $grpid = shift();
	my $gbFs = $GENBLASTDIR . $GENBLASTSUBDIR . $grpid . "*";
	`rm -f $gbFs`;
	return;
}

#+ -------------- for GFF3 gene model comparison
sub get_existing_gms {
# called my $ensgdidR =  Genblast::get_existing_gms($gffppF); 
	#reads Ensembl preprocessed gff file 
	my $gffppF = shift();
	die "gff3 pp file $gffppF not there\n" unless (-e $gffppF); 
	my %ensgd_id;
	open (IN, $gffppF);	
	while (<IN>) {
		chomp;
		my ($gd, $id) = split("\t", $_);
		$ensgd_id{$gd} = $id;
	}
	close IN;
	return \%ensgd_id;
}

sub genedet_genblast {
# called my $gbgd =  Genblast::genedet_genblast($gbgffF); 
	#reads genblast gff file and returns gene details text
	my $gbgffF = shift();
	die "Genblast gff file $gbgffF not there\n" unless (-e $gbgffF); 
	my ($gbgsR, $gbgeneR) = parse_gff3("genblast", $gbgffF);
	my $gbg = $$gbgsR[0];
	my $gbgd = gene_details_txt($gbgeneR->{$gbg});
	return $gbgd;
}

sub parse_gff3 {
	my $source = shift(); 
	my $gffF = shift(); 
	my $gbcds;
	my $uselimit = ($gbcds = shift());
	my $gffio = Bio::Tools::GFF->new(-file => $gffF, -gff_version => 3);
	my %gene;
	my @gs;
	# loop over the input stream
	while (my $feat = $gffio->next_feature()) {
		if ($feat->primary_tag eq $GFFPRIMARYTAGS{$source}{cds}) {
			next if ($uselimit && ! $gbcds->overlaps($feat, 'strong'));
			die "No ID tag in $gffF transcript?\n" unless $feat->has_tag ('ID');
			my @a = $feat->get_tag_values('ID'); my $g = $a[0];
			if (exists $gene{$g}) {
				warn "gene $g already exists. $gffF transcript\n" ;
			} else {
				push @gs, $g;
			}
			$gene{$g}{cds} = $feat;
		} elsif ($feat->primary_tag eq $GFFPRIMARYTAGS{$source}{"exon"}) { 
			die "No Parent tag in $gffF exon?\n" unless $feat->has_tag ('Parent');
			my @a = $feat->get_tag_values('Parent'); my $g = $a[0];
			unless (exists $gene{$g}) {
				next if ($uselimit);
				warn "gene $g does not exist. $gffF $gffF \n";
				push @gs, $g;
			}
			push @{ $gene{$g}{exons}}, $feat;	  
		}
	}
	$gffio->close();
	return (\@gs, \%gene);
}

sub gene_details_txt {
	# returns a scalar txt which completely defines gene model cds
	my $g = shift();
	my @dets = ($g->{cds}->seq_id(), $g->{cds}->strand());
	my @exons = @{ $g->{exons}};
	push @dets, scalar @exons;
	my @exonstxts;
	foreach my $exon (@exons) {
		push @exonstxts, $exon->start(). "..". $exon->end();
	}
	push @dets, join(" ", @exonstxts);
	return join(";", @dets);
}
#- -------------

1;
