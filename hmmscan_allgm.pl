#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;
use HMMs;

# run fasta peptide file of all sp gene models against current HMM database of all grps

# defaults
my $DEF_DBOPT = "final_hmm";
my $DEF_SSOPT = "all";

## get inputs
my $qsp = $ARGV[0] || die "first argument must be query sp or \'none\' to make db\n";
my $dbopt = $ARGV[1] || $DEF_DBOPT; # "HMM1"|"final_hmm";
my $ssopt = $ARGV[2] || $DEF_SSOPT; # all or subset file e.g. "genes_with_evidence.txt";
die "second argument must be must be \'HMM1\' or \'final_hmm\' not $dbopt\n" unless ($dbopt eq "HMM1" || $dbopt eq "final_hmm");
# per sp status file (GM steps ; must be all complete). defines @sps
my ($gmfinished, $spstatR) = Pipeline::read_gmsteps_status();	
die "All genemodel need to be complete to run this\n" unless ($gmfinished);

my $makedbonly = ($ARGV[0] eq "none");
my $qF; 
unless ($makedbonly) {
	my $q_is_grass = 0;
	foreach my $gsp (@sps) {
		$q_is_grass = ($qsp eq $gsp);
		last if ($q_is_grass);
	}
	die "query must be grass in HMM1 mode\n" unless ($q_is_grass || $dbopt eq "final_hmm");
	die "No subset allowed in HMM1 mode\n" unless ($ssopt eq "all" || $dbopt eq "final_hmm");
	if ($q_is_grass) {
		print "query sp $qsp is a grass. HMM mode: $dbopt group set: $ssopt\n";
		$qF = Pipeline::gmfasta_filename( $qsp );
	} else {
		print "query sp $qsp is a nongrass. HMM mode: $dbopt group set: $ssopt\n";
		my ($qspsR, $ngsp_infoR) = Pipeline::nongrass_info();
		my %ngsp_info = %$ngsp_infoR;
		$qF = $ngsp_info{$qsp}{pep_fasta};
	}
	die "query file $qF not there\n" unless (-e $qF );
}

my ($grpid_statusR, $fsstatF) = Pipeline::get_fsteps_status("all"); # all for all and subset
# make a HMMER db of all groups with current hmm
# called my ($db, $dbgrpidsR) = HMMs::make_db($grpid_statusR, "HMM1"|"final_hmm", "all"); 
# or
# called my ($db, $dbgrpidsR) = HMMs::make_db($grpid_statusR, "final_hmm", $ssF); 
my ($db, $dbgrpidsR) = HMMs::make_db($grpid_statusR, $dbopt, $ssopt);
my $nprofile=@$dbgrpidsR;
print "\nThere are $nprofile HMMER profiles in $db\n";
exit(0) if ($makedbonly);
#-

#+--------------  run hmmscan 
my $hmmFsR = HMMs::run_query_hmmscan($qsp, $qF, $db);
my $hitsR = HMMs::get_hmm_scores($hmmFsR);
print "$qsp hmmscan complete\n";
#-

my $outputmsg;
if ($dbopt eq "HMM1") {
	$outputmsg = HMMs::output_hmmscan_hits($hitsR, $grpid_statusR, $qsp, $db, "topq");
} else {
	$outputmsg = HMMs::output_hmmscan_hits($hitsR, $grpid_statusR, $qsp, $db, "all");
}
print $outputmsg;

exit(0);

	
		
	

