#!/usr/bin/env perl

##################################################
# universal_grass_peps pipeline v1.3             #
# Rowan Mitchell 2023                            #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################


use strict;
use warnings;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;
use Ortho;

# constants
#+ for reading in Ensembl ortholog tables (manually downloaded)
# inputs ortholog tables and output essential info to simpler file for each ref spp
my @HCONTS = ("Gene stable ID", "Transcript stable ID", "gene stable ID", "protein or transcript stable ID", "Query protein or transcript ID", 
"Last common ancestor with ", "homology type", 	"%id. target", "%id. query", "orthology confidence"  );		
# full set is ("Gene stable ID", "Transcript stable ID", "gene stable ID", "protein or transcript stable ID", "Query protein or transcript ID", 
# "Last common ancestor with ", "homology type", 	"%id. target", "%id. query", 	"Whole-genome alignment coverage", 	
# "orthology confidence"  );		
#-


my ($gmfinished, $sps_ref) = Pipeline::read_gmsteps_status();	


#+ get redundant id mapping to nr id used throughout	
my %red_id;
open (RED, "$REDSEQF"); # get redundant ids
while (<RED>) {
	chomp;
	my ($sp, $redid, $id) = split("\t", $_);
	$red_id{$redid} = $id;
}
close RED;
#- 
my %refsp_clust_peps; # clustering of peps into groups from self blastp
my %refsp_pep_clust; # clustering of peps into groups from self blastp

foreach my $refsp (@REFSPS) {
	print "\nRef species $refsp\n";
	my $outorthoF = $REFSP_SHTN{$refsp} . $OUTORTHOEXT;
	($refsp_clust_peps{$refsp}, $refsp_pep_clust{$refsp}) = Ortho::get_clusters($refsp); # all member peps as hash clust -> [@members]

	my %clust_orth; # reorganised orths for clusters ; main output
	
	foreach my $sp (@sps) {
		next if ($sp eq $refsp);
		print "\tspecies $sp\n";
		my $orthoF = $REFSP_SHTN{$refsp} . "-" . $sp . ".txt";
		die "File $orthoF not found\n" unless (-e $orthoF);
		open (IN, $orthoF);
		#+ set input cols from header
		my $htxt = <IN>;
		chomp $htxt;
		my @hs = split("\t", $htxt);
		my %hcont_i;
		for (my $i=0; $i<@hs; $i++) {
			foreach my $hcont (@HCONTS) {
				 if ($hs[$i] =~ /$hcont/) {$hcont_i{$hcont} = $i; last;}
			}
		}
		foreach my $hcont (@HCONTS) {
			die "$hcont not found in header $htxt from $orthoF\n" unless (exists $hcont_i{$hcont});
		}
		#-
	my %allorth;
	LINE: while (<IN>) {
			chomp;
			my @a = split("\t", $_);	
			for (my $i=0; $i<@hs; $i++) {
				next LINE if (!defined $a[$i] || length $a[$i] == 0);
			}
			my $q = $a[ $hcont_i{"Query protein or transcript ID"}]; # $q is query pep id from refsp
			$q = $red_id{ $q} if (exists $red_id{ $q}); # set to nr id
			my $target = $a[ $hcont_i{"protein or transcript stable ID"}];
			$target = $red_id{ $target} if (exists $red_id{ $target}); # set to nr id
			my %rec;
			$rec{t_pcid} = $a[ $hcont_i{"%id. target"}]; 
			$rec{q_pcid} = $a[ $hcont_i{"%id. query"}]; 
			$rec{confidence} = $a[ $hcont_i{"orthology confidence"}]; 
			$allorth{$q}{$target} = \%rec;
		}
		close IN;
		foreach my $q (keys %allorth) {
			my $clust = $refsp_pep_clust{$refsp}{$q}; 
			my @peps = @{ $refsp_clust_peps{$refsp}{$clust}};
			next unless ($q eq $peps[0] ); # only get orthos for top ranked (longest) pep in cluster 
			foreach my $t (keys %{ $allorth{$q}} ) { 
				$clust_orth{$clust}{$sp}{$t} = $allorth{$q}{$t};
			}
		}
		
	} # end each sp
	
	open (OUT, ">$outorthoF");
	my @h = ("cluster_id", @sps);
	print OUT join("\t", @h) . "\n";
	my $nout=0;
	foreach my $clust (sort keys %clust_orth) {
		my @row = ($clust);
		foreach my $sp (@sps) {
			if (exists $clust_orth{$clust}{$sp}) {
				my $peps_ref = Ortho::sort_orths( $clust_orth{$clust}{$sp});
				my @cs;
				foreach my $pep (@$peps_ref) {
					my %h = %{ $clust_orth{$clust}{$sp}{$pep} };
					push @cs, join(" ", ($pep, @h{@ORTHKEYS}));
				}
				push @row, join(";", @cs);
			} else {
				push @row, "";
			}
		}
		print OUT join("\t", @row) . "\n";
		$nout++;
	}
	close OUT;
	print "Wrote all orths of $nout grps for $refsp to $outorthoF\n";
			
}

exit(0);

