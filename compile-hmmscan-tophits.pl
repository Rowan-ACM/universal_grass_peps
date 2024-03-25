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

# reads in all (available) nongrass tophits from hmmscans against final_db, compiles into 1 file and gives output stats
# adds specifcity info columns to final db table

#+ constants
my $DB = "final_db";
my $CHTF = "compiled-nongrass-hmmscan-tophits.txt";
my $FDBWSF = "final_db_with_spec.txt";
## define hierarchical categories : 
# "non-grass" is all hits
# "non-commelinid" is hits excluding commelinid monocot
# "non-monocot" is hits excluding any monocot
my @MCCATS = ("non-grass", "non-commelinid", "non-monocot");
# corresponding specificity categories
my %MC_SPCAT = ("non-grass"=>"grass specific", "non-commelinid"=>"commelinid specific", "non-monocot"=>"moncot specific");
my $NMC = 2;
my $DEFCUTOFFS = 0.25; # default cut-off define specific
#-

#+ following is just to get set of grpids used to make db and check that this fits with those found in hmmscan results
my $gmfinished = Pipeline::read_gmsteps_status();	
die "All genemodel need to be complete to run this\n" unless ($gmfinished);
my ($grpid_statusR, $statF) = Pipeline::get_fsteps_status("all"); # get final_db status file
my @dbgrpids = sort keys %$grpid_statusR;
my $nprofile = @dbgrpids;
print "\nThere are $nprofile HMMER profiles in $DB\n";
#-

my ($qspsR, $ngsp_infoR) = Pipeline::nongrass_info();
my %ngsp_info = %$ngsp_infoR;

my %table; # all hmmscan rel score values to be output
my $nspdone=0;
foreach my $qsp (@$qspsR) {
	my $hmmscanthF = HMMs::get_hmmscan_hits_filename($qsp, $DB, "all");
	unless (-e $hmmscanthF) {
		warn "\tHit file not found for $qsp\n\t$hmmscanthF\n";
		next;
	}
	open (IN, $hmmscanthF) || die "cannot open $hmmscanthF\n";
	<IN>;
	while (<IN>) {
		chomp;
	# cols: query   grpid   hit_score       maxposs_score   hit_relscore    lowestmem_relscore      diff_lowestmem_hit_relscore(= specificity)
		my ($q, $grpid, $ngscore, $maxscore, $ngrelscore, $lgrelscore, $difrelscore) = split("\t", $_);
		die "Unrecognised group id $grpid in $hmmscanthF\n" unless (exists $grpid_statusR->{$grpid});
		#+ get top hit (lowest specificity) for this grpid, sp
		$table{$grpid}{$qsp} = $difrelscore unless (exists $table{$grpid}{$qsp} && $table{$grpid}{$qsp} < $difrelscore); 
		#-
	}
	close IN;
	$nspdone++;
}
my $nqsp = @$qspsR;
print "\n\nFound output for $nspdone spp out of $nqsp\n";

#+ ------------------------ output
my %mindiff; 
open (OUT, ">$CHTF");
open (FWS, ">$FDBWSF"); # final_db with added spec
my @taxa;
foreach my $qsp (@$qspsR) {
	push @taxa, $ngsp_info{$qsp}{taxon};
}
print OUT join("\t", ("hmmscan diff rel scores (lowest grass - highest nongrass)", @taxa)) . "\n"; # header info #1
my @mincts; # for OUT
my @spechs; # specificity col headers for FWS
foreach my $m (@MCCATS) {
	push @mincts, "min $m"; 
	push @spechs, "S $m";
	push @spechs, "sp of top hit $m";
}
push @spechs, "categorisation using S<$DEFCUTOFFS";
print OUT join("\t", ("group_id", @$qspsR, @mincts)) . "\n"; # column titles
print FWS join("\t", ("group_id", @spechs, @GRPDETSUM, @sps)) . "\n"; # FWS column titles
my %grpid_status = %$grpid_statusR;

foreach my $grpid (@dbgrpids) {
	my @row = ($grpid);		
	my @fwsrow = ($grpid);
	foreach my $qsp (@$qspsR) {
# get diff values		
		if (exists $table{$grpid}{$qsp}) {
			push @row, $table{$grpid}{$qsp};
			my $ncat; # how many categories to check against
			if ($ngsp_info{$qsp}{taxon} =~ /Commelinid/i) {
				$ncat=1;
			} elsif ($ngsp_info{$qsp}{taxon} =~ /Monocot/i) {
				$ncat=2;
			} else {
				$ncat=3;
			}
			for (my $i=0;$i<$ncat;$i++) {
				my $mccat = $MCCATS[$i];
				$mindiff{$mccat}{$grpid} = {"diff"=>$table{$grpid}{$qsp}, "species"=>$qsp} if 
				  (! exists $mindiff{$mccat}{$grpid} || $table{$grpid}{$qsp} < $mindiff{$mccat}{$grpid}{"diff"});
			}
		} else {
			push @row, "";
		}
	}
	my $speccat;
	foreach my $mccat (@MCCATS) {
		unless (exists $mindiff{$mccat}{$grpid}) { # no hit so min_diff (S) is equal to grp "mingrp_relscore"
			$mindiff{$mccat}{$grpid}{"diff"} = $grpid_statusR->{$grpid}{"current_hmm_details"}{"mingrp_relscore"};
			$mindiff{$mccat}{$grpid}{"species"} = "no_hit";
		}
		my $cell = sprintf "%.2f %s%s%s", $mindiff{$mccat}{$grpid}{"diff"}, "(", $mindiff{$mccat}{$grpid}{"species"}, ")";
		push @row, $cell;
		push @fwsrow, $mindiff{$mccat}{$grpid}{"diff"};
		push @fwsrow, $mindiff{$mccat}{$grpid}{"species"};
		$speccat = $MC_SPCAT{$mccat} if ($mindiff{$mccat}{$grpid}{"diff"} > $DEFCUTOFFS && ! defined $speccat);
	}
	$speccat = "non-specific" unless (defined $speccat);
	print OUT join("\t", @row) . "\n";
	push @fwsrow, $speccat;
	my @gds = @{ Pipeline::_group_details_texts($grpid_status{$grpid}{"current_hmm_details"})};
	@fwsrow = (@fwsrow, @gds); # append row with grp hmm details
	print FWS join("\t", @fwsrow). "\n";
}
close FWS;
close OUT;
print "wrote non-grass hit output to $CHTF and added specificity to final_db universal set in $FDBWSF\n";


#--- specificity stats
print "\nSummary stats for the total $nprofile groups:\n";
my %mostspec;
foreach my $sc (@MC_SPCAT{@MCCATS}) {
	$mostspec{$sc} = {"diff"=>0, "grpid"=>"none", "species"=>"none"};
}
foreach my $grpid (@dbgrpids) {
	foreach my $mccat (@MCCATS) {
		my $sc = $MC_SPCAT{$mccat};
		if ($mindiff{$mccat}{$grpid}{"diff"} < 1) {
			my $diff = $mindiff{$mccat}{$grpid}{"diff"};
			$mostspec{$sc} = {"diff"=>$diff, "grpid"=>$grpid} if ($diff > $mostspec{$sc}{"diff"});
		}
	}
}
print join("\t", ("most spec grpid", @MC_SPCAT{@MCCATS})) . "\n";
my @row = ("grpid:");
foreach my $sc (@MC_SPCAT{@MCCATS}) {
	my $txt = sprintf "%s (%.2f)", $mostspec{$sc}{"grpid"}, $mostspec{$sc}{"diff"};
	push @row, $txt;
}
print join("\t", @row) . "\n\n";
# 2. find most specific grpidscount numbers passing different upper limits
for (my $limit=0.0; $limit<0.5; $limit+=0.1) {
	my %nspec;
	foreach my $grpid (@dbgrpids) {
		foreach my $mccat (@MCCATS) {
			my $sc = $MC_SPCAT{$mccat};
			$nspec{$sc}++ if ($mindiff{$mccat}{$grpid}{"diff"}	> $limit);
		}
	}
	my @statrow = ("> $limit\:");
	foreach my $sc (@MC_SPCAT{@MCCATS}) {
		next unless (exists $nspec{$sc});
		my $pc = 100* $nspec{$sc} / $nprofile; my $pctxt = sprintf "%.0f", $pc;
		push @statrow, "$nspec{$sc} ($pctxt\%)";
	}
	print join("\t", @statrow) . "\n";
}


