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

# reads in grasshits from hmmscan files for all grass sp
# then finds all associated peptides 
# outputs info for all grass sp as a table of counts for :
# genblast gene models (does not use hmmscan data) ngbgm, ngrp_ass_peps, max_ass_peps, tot_ass_peps

# also output info on supergrps

# constants
#+ main output files to summarise final database
my $ASSOUTF = "asspeps.txt";
my $ABMOUTF = "ass_bettermatches.txt"; # for ass peps that would apparently be better matches than existing member. for info only
my $SUPERGRPOUTF = "supergrps.txt";
my $SPTABF = "allgrass_sp_table.txt";
my @SPTABCOLS = ("ngbgm", "ngrp_ass_peps", "max_ass_peps", "tot_ass_peps");
my %SPTABCOL_H = (ngbgm=>"num genBlastG gene models", ngrp_ass_peps=>"num grps with ass peps"
                  , max_ass_peps=>"max num ass peps in 1 grp", tot_ass_peps=>"total num ass peps");
#+-

# per sp status file (GM steps ; must be all complete). defines @sps
my $gmfinished = Pipeline::read_gmsteps_status();	
die "All genemodel need to be complete to run this\n" unless ($gmfinished);
my ($grpid_statusR, $fsstatF) = Pipeline::get_fsteps_status("all"); 

my ($db, $dbgrpidsR, $sp_grpid_asspepR, $sp_grpid_bettermatchesR, $grpid_qgrpid_spR) = HMMs::get_nonmem_hitscores($grpid_statusR, "final_hmm", "all", $MINRELSC);
my @dbgrpids = sort keys %$dbgrpidsR;
my $ngrp = @dbgrpids;
print "$ngrp grpids used in db $db\n";
#+ output files of hits for each grp and sp
my $nrow;
$nrow = HMMs::output_nonmem_scores($ASSOUTF, $sp_grpid_asspepR, $dbgrpidsR);
print "Wrote $nrow rows associated peps to file $ASSOUTF\n";
$nrow = HMMs::output_nonmem_scores($ABMOUTF, $sp_grpid_bettermatchesR, $dbgrpidsR);
print "Wrote $nrow rows better matches to file $ABMOUTF\n";

my %sp_table; # all per sp summary info for output
#+ summarise genblast info
my $ngbtot=0;
my %grpid_nogb; @grpid_nogb{@dbgrpids}=(1) x @dbgrpids; 
print "\n------ Genblast info -----------\nsp\tn genblast gene models\n";
my %sp_grpid_isgblast;
foreach my $grpid (@dbgrpids) {
	foreach my $sp (@sps) {
		my $member = $grpid_statusR->{$grpid}{"current_hmm_details"}{"members"}{$sp};
		my $membid = $member->{"id"};	
		$sp_grpid_isgblast{$sp}{$grpid} = $membid =~ /^genblast/;
	}
}
foreach my $sp (@sps) {
	my $ngb = 0;
	foreach my $grpid (@dbgrpids) {
		if  ( $sp_grpid_isgblast{$sp}{$grpid} ) {
			$ngb++;
			$grpid_nogb{$grpid} = 0;
		}
	}
	print "$sp\t$ngb\n";
	$sp_table{$sp}{"ngbgm"} = $ngb;
	$ngbtot += $ngb;
}	
print "total\t$ngbtot\n";
my $ngrpnogb = 0; foreach my $grpid (@dbgrpids) {$ngrpnogb++ if ($grpid_nogb{$grpid}); }
print "\n$ngrpnogb out of $ngrp grpids have no genblast members\n";
print "--------------------------------\n\n";
#-

#+ generate summary info
foreach my $qsp (@sps) {

	# summary stats on ass pep
	my $nmultiass=0;
	my $ntotass=0;
	my $maxnass=0;
	my $nbm=0;
	my $nnf_nobm=0;
#	my @notfound;
	foreach my $grpid (@dbgrpids) {	
		my $nass = 0; 
		if (exists $sp_grpid_asspepR->{$qsp}{$grpid}) {
			$nass = scalar keys %{ $sp_grpid_asspepR->{$qsp}{$grpid}}; # number of non-member queries associated to this grp 
		}
		$nmultiass++ if ($nass > 0);
		$maxnass = $nass if ($nass > $maxnass);
		$ntotass += $nass;
		if (exists $sp_grpid_bettermatchesR->{$qsp}{$grpid}) {
			$nbm++; # number of non-member better matches than grp member  
		# } else {
			# $nnf_nobm++ unless ($grpid_found{$grpid} || $sp_grpid_isgblast{$qsp}{$grpid}); # member not found but no better match?
		}	
#		push @notfound, $grpid unless ($grpid_found{$grpid} || $sp_grpid_isgblast{$qsp}{$grpid});
		
	}

	print "\nSummary of associated peps for species **$qsp**:\n"; 
	print "\t$nmultiass have at least 1 associated peptide\n";
	$sp_table{$qsp}{ngrp_ass_peps} = $nmultiass;
	print "\t$maxnass is max number associated peptides to grp  \n";
	$sp_table{$qsp}{max_ass_peps} = $maxnass;
	my $npep = $ngrp + $ntotass;
	print "\t$ntotass total associated peptides so $npep peptides are in or associated to grps \n"; 
	$sp_table{$qsp}{tot_ass_peps} = $ntotass;
	# my $nnf = @notfound;
	# if ($nnf == 0) {
		# print "\nMember peps found as expected for all groups\n";
	# } else {
		# print "\nNumber groups with: missing members $nnf, better matches $nbm, missing member no better match $nnf_nobm\n";
	# }

}
#- 

#output sp summary table
open (OUT, ">$SPTABF");
my @hs = ("grass sp", @SPTABCOL_H{@SPTABCOLS});
print OUT join("\t", @hs) . "\n";
foreach my $sp (@sps) {
	my @row =($sp);
	foreach my $col (@SPTABCOLS) {
		push @row, $sp_table{$sp}{$col};
	}
	print OUT join("\t", @row) . "\n";
}
close OUT;

print "\nDefining supergroups:\n";
my %grpid_qgrpid_sp = %{$grpid_qgrpid_spR};
#+ supergroups definition
my %supergrp_grpid; # 1 to many 
my %grpid_supergrp; # 1 to 1
foreach my $grpid (sort keys %grpid_qgrpid_sp) {	
	#+ all matching other grpids 
	

QGRPID: foreach my $qgrpid (sort keys %{ $grpid_qgrpid_sp{$grpid} } ) {	
		
		
		foreach my $sp (@sps) {
			next QGRPID unless (exists $grpid_qgrpid_sp{$grpid}{$qgrpid}{$sp} || $sp_grpid_isgblast{$sp}{$grpid}); 
		}
					
		# this grpid - qgrpid pair have all members matching
		if (exists $grpid_supergrp{ $grpid }) {
			if (exists $grpid_supergrp{ $qgrpid }) { # both exist
				next QGRPID if ($grpid_supergrp{ $grpid } == $grpid_supergrp{ $qgrpid }); # already consistent so finish
				# inconsistent- means another grp matching grpid did not fully match qgrpid which matches different other 
				# because grpid - qgrpid check is stringent, this means supergrps are v close so combine
				#+ combine supergrps
				my ($keepgrp, $delgrp);
				if ($grpid_supergrp{ $grpid } < $grpid_supergrp{ $qgrpid }) { # keep lower number ID (arbitrary)
					$keepgrp = $grpid;
					$delgrp = $qgrpid;
				} else {
					$keepgrp = $qgrpid;
					$delgrp = $grpid;
				}
				my @mvgrps = keys %{ $supergrp_grpid{ $grpid_supergrp{ $delgrp }}};
				delete $supergrp_grpid{ $grpid_supergrp{ $delgrp }};
				foreach my $mvgrp (@mvgrps ) {
					$supergrp_grpid{ $grpid_supergrp{ $keepgrp }}{$mvgrp} = 1;
					$grpid_supergrp{ $mvgrp } = $grpid_supergrp{ $keepgrp };
				}							
				#-  combine supergrps
			} else {
				#  add qgrpid to supergrp that grpid is in
				$supergrp_grpid{ $grpid_supergrp{ $grpid }}{$qgrpid} = 1;
				$grpid_supergrp{ $qgrpid } = $grpid_supergrp{ $grpid };
			}
		} elsif (exists $grpid_supergrp{ $qgrpid }) { 
			#  add grpid to supergrp that qgrpid is in
			$supergrp_grpid{ $grpid_supergrp{ $qgrpid }}{$grpid} = 1;
			$grpid_supergrp{ $grpid } = $grpid_supergrp{ $qgrpid };
		} else { # neither exist so create new supergrp 
			my @supergpids = sort {$a <=> $b} keys %supergrp_grpid;
			my $supergrpid;
			if (@supergpids eq 0) {
				$supergrpid = 1;
			} else {
				$supergrpid = $supergpids[$#supergpids] + 1; # new $supergrpid is 1 more than last element of sorted keys 
			}
			$supergrp_grpid{ $supergrpid }{$grpid} = 1;
			$supergrp_grpid{ $supergrpid }{$qgrpid} = 1;
			$grpid_supergrp{ $grpid } = $supergrpid;
			$grpid_supergrp{ $qgrpid } = $supergrpid;
		}
	}
}
#- supergroups definition


# get ids sorted by descending numbers of grps
my @supergpids = sort {scalar keys %{$supergrp_grpid{$b}} <=> scalar keys %{$supergrp_grpid{$a}}} keys %supergrp_grpid;
open (OUT, ">$SUPERGRPOUTF");
print OUT "supergrpID\tn groups\tgrpids\n";
foreach my $supergrpid (@supergpids) {
	my @grpids = keys %{$supergrp_grpid{ $supergrpid }}; 
	my $ng = scalar @grpids;
	print OUT "\#$supergrpid\t$ng\t@grpids\n";
}
close OUT;
my $nsupergrp = scalar keys %supergrp_grpid;
print "\t$nsupergrp supergroups identified. Written to $SUPERGRPOUTF\n";
exit (0);

