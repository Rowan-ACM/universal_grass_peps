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

# run fasta peptide file of all sp gene models against current HMM database of all grps

# defaults
my $DEF_CHMM = "final_hmm";
my $DEF_SSOPT = "all";

## get inputs
my $qsp = $ARGV[0] || die "first argument must be query sp or \'none\' to make db\n";
my $chmm = $ARGV[1] || $DEF_CHMM; # "HMM1"|"final_hmm";
my $ssopt = $ARGV[2] || $DEF_SSOPT; # all or subset file e.g. "genes_with_evidence.txt";
die "second argument must be must be \'HMM1\' or \'final_hmm\' not $chmm\n" unless ($chmm eq "HMM1" || $chmm eq "final_hmm");
# per sp status file (GM steps ; must be all complete). defines @sps
my $gmfinished = Pipeline::read_gmsteps_status();	
die "All genemodel need to be complete to run this\n" unless ($gmfinished);

my ($grpid_statusR, $statF) = Pipeline::get_fsteps_status("all"); # get current status file
my ($db, $dbgrpidsR) = HMMs::dbname_grpset($grpid_statusR, $chmm, $ssopt);

my $makedbonly = ($ARGV[0] eq "none");
my $qF; 
if ($makedbonly) {
	my $dbmsg = HMMs::make_db($db, $dbgrpidsR);
	print "\n$dbmsg\n";
	exit(0);
} else { # set qF
	my $q_is_grass = 0;
	foreach my $gsp (@sps) {
		$q_is_grass = ($qsp eq $gsp);
		last if ($q_is_grass);
	}
	die "query must be grass in HMM1 mode\n" unless ($q_is_grass || $chmm eq "final_hmm");
	die "No subset allowed in HMM1 mode\n" unless ($ssopt eq "all" || $chmm eq "final_hmm");
	if ($q_is_grass) {
		print "query sp $qsp is a grass. HMM mode: $chmm group set: $ssopt\n";
		$qF = Pipeline::gmfasta_filename( $qsp );
	} else {
		print "query sp $qsp is a nongrass. HMM mode: $chmm group set: $ssopt\n";
		my ($qspsR, $ngsp_infoR) = Pipeline::nongrass_info();
		$qF = $ngsp_infoR->{$qsp}{pep_fasta};
	}
	die "query file $qF not there\n" unless (-e $qF );
}
#-

#+--------------  run hmmscan 
my $hmmFsR = HMMs::run_query_hmmscan($qsp, $qF, $db);
my $hitsR = HMMs::get_hmm_scores($hmmFsR);
print "$qsp hmmscan complete\n";
#-

my $outputmsg;
if ($chmm eq "HMM1") {
	$outputmsg = HMMs::output_hmmscan_hits($hitsR, $grpid_statusR, $qsp, $db, "topq");
} else {
	$outputmsg = HMMs::output_hmmscan_hits($hitsR, $grpid_statusR, $qsp, $db, "all");
}
print $outputmsg;
exit(0);

	
		
	

