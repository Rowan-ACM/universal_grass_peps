#!/usr/bin/env perl

##################################################
# universal_grass_peps pipeline v1.4             #
# Rowan Mitchell 2020-2024                       #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################

# * (no change from v1.3)

use strict;
use warnings;
use Bio::SeqIO;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;
use Ortho;

# writes group members identified in prev steps to fasta files for all groups that pass cut-off for min number of grass spp
# also writes alt poss seqs from gene models non-assigned to groups for later seaching by hmmscan

my $gmfinished = Pipeline::read_gmsteps_status();	

my @spstodo;
foreach my $sp (@sps) {
	push @spstodo, $sp unless ($sp_stat{$sp}{"gmstep_completed"} eq $COMPLETETXT);
}
my $ntodo = @spstodo;

# get all orthologs and species 
my ($OTF) = Ortho::read_orth_table();
print "Write fasta files for orthos in $OTF. $ntodo spp to do\n\n";

my $firstsp=1; # used to check if file is already there
foreach my $sp (@spstodo) {
	print "\nWriting ortho table entries for $sp\n";
	my $gmfaF = Pipeline::gmfasta_filename( $sp );
	my ($idsR, $id_seqobjR) = Pipeline::read_fasta( $gmfaF );
	my %id_seqobj = %{ $id_seqobjR};

	my $ngrp=0;
	foreach my $grpid (sort keys %grp_table) {
		if (exists $grp_table{$grpid}{$sp}) {
			my ($id, @altids) = @{ $grp_table{$grpid}{$sp} };
			die "seq obj does not exist for $sp : $id \n" unless (exists $id_seqobj{$id});
			my $faFsR = Pipeline::group_fasta_fnames($grpid); # get fasta filenames 
			my $hmm0faF = $faFsR->{"HMM0"};	
			my $altfaF = $faFsR->{"HMM0alt"};	
			die "File already exists $hmm0faF - delete and re-run\n" if ($firstsp && -e $hmm0faF);
			die "File already exists $altfaF - delete and re-run\n" if ($firstsp && -e $altfaF);
			Pipeline::add_seq_to_grpfasta($sp, $hmm0faF, $id_seqobj{$id} );
			my %altid_seqobj;
			@altid_seqobj{ @altids } = @id_seqobj{@altids};
			my $nseq = Pipeline::append_fasta($altfaF, \%altid_seqobj );
			$ngrp++;
		}	
	}
	$sp_stat{$sp}{"gmstep_completed"} = $COMPLETETXT;
	Pipeline::update_gmsteps_status();
	print "\twrote $sp seqs to fasta files for $ngrp groups\n";
	$firstsp=0;
}

# generate empty per-grp status file to be used by subsequent steps in pipeline
Pipeline::create_fsteps_status_file();
print "\nEmpty file of status for each group written to $FSSTATFALL\n";
exit(0);

