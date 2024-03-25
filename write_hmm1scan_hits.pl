#!/usr/bin/env perl

##################################################
# universal_grass_peps pipeline v1.4             #
# Rowan Mitchell 2020-2024                       #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################


use strict;
use warnings;
use Bio::SeqIO;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;
use HMMs;

# reads in grasshits from hmmscan files for all grass sp at HMM1 stage 
# outputs info as a table and write fasta files for each grp
# similar to EnsemblPlants/orthologs/write-all-grp-fa.pl
my $RELSCCUTOFF = 0.5; # saves time by eliminating hits that are too low score to be plausible alternatives

my $gmfinished = Pipeline::read_gmsteps_status();	# per sp status file (GM steps ; must be all complete). defines @sps
die "All genemodel need to be complete to run this\n" unless ($gmfinished);
my ($grpid_statusR, $fsstatF) = Pipeline::get_fsteps_status("all"); 

# callled my ($db, $dbgrpidsR, $sp_grpid_asspepR, $sp_grpid_bettermatchesR, $grpid_qgrpid_spR) = HMMs::get_nonmem_hitscores($grpid_statusR, $chmm, $hitsopt, $cutoff);
my ($db, $dbgrpidsR, $sp_grpid_asspepR, $sp_grpid_bettermatchesR, $grpid_qgrpid_spR) = HMMs::get_nonmem_hitscores($grpid_statusR, "HMM1", "topq", $RELSCCUTOFF);
my @dbgrpids = sort keys %$dbgrpidsR;
my $ngrp = @dbgrpids;
print "$ngrp grpids used in db $db\n";

#+ output file of hits for each grp and sp
my $nrow = HMMs::output_nonmem_scores($HMM1HITSF, $sp_grpid_bettermatchesR, $dbgrpidsR);
print "Wrote $nrow rows with at least 1 better match to $HMM1HITSF\n";

my $firstsp=1; # used to check if file is already there
foreach my $sp (@sps) {
	print "\nGetting all fasta for $sp\n";
	my $gmfaF = Pipeline::gmfasta_filename( $sp );
	my ($idsR, $id_seqobjR) = Pipeline::read_fasta( $gmfaF );
	my %id_seqobj = %{ $id_seqobjR};

	my $ngrp=0;
	foreach my $grpid (@dbgrpids) {
		if (exists $sp_grpid_bettermatchesR->{$sp}{$grpid}) {
			my $faFsR = Pipeline::group_fasta_fnames($grpid); # get fasta filenames 
			my $hmm1altfaF = $faFsR->{"HMM1alt"};	
			die "File already exists $hmm1altfaF  - delete and re-run\n" if ($firstsp && -e $hmm1altfaF);
			my @altids = keys %{ $sp_grpid_bettermatchesR->{$sp}{$grpid} };	
			my %altid_seqobj;
			@altid_seqobj{ @altids } = @id_seqobj{@altids};
			my $nseq = Pipeline::append_fasta($hmm1altfaF, \%altid_seqobj );
			$ngrp++;
		}	
	}
	print "\twrote $sp seqs of hits to fasta files for $ngrp groups\n";
	$firstsp=0;
}

exit (0);

