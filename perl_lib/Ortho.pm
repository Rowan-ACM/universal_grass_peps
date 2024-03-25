##################################################
# universal_grass_peps pipeline v1.4             #
# Rowan Mitchell 2020-2024                       #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################



package Ortho;
use strict;
use warnings;
use Pipeline;

#
# following are passed to subroutines as global variables (large or unchanged). Exported from Orth
# our %refsp_clust_peps; # clustering of peps into groups from self blastp clust -> [@peps]
# our %refsp_pep_clust; # clustering of peps into groups from self blastp	pep -> clust
# our %refsp_clust_allorth; # all orthologs for each refsp clust
# our %pepset; # all seed peps fron ref spp
# our %orth_topclust; # top matching clust for each orth
# our %grp_table; # Table of orths that is main output of script

sub get_clusters {
# read in clusters of refsp peps, clustered by self blastp
	my $refsp = shift();
	my $clustF = $BLASTPDIR . $refsp . $CLUSTFEXT;
	die "$refsp cluster file $clustF not found\n" unless (-e $clustF);
	my %clust_peps; # cluster name  -> array members
	my %pep_clust; # member -> cluster name
	open (CF, "$clustF");
	<CF>;
	while (<CF>) {
		chomp;
		my ($clustname, @peps) = split("\t", $_);
		my ($b, $clust) = split(":", $clustname);
		$clust_peps{$clust} = [@peps];	
		foreach my $pep (@peps) {$pep_clust{$pep} = $clust;}
	}
	close CF;
	return (\%clust_peps, \%pep_clust);
}

sub read_in_Zm_blastp_Os {
	# usage my $osclust_blastp_zmclustR = Ortho::read_in_Zm_blastp_Os();
	#+ read Zm blastp hits for Os and use this to provide 1:1 mapping of os cluster -> zm cluster 	
	my %osclustcheck;
	my %zmclust_tophit;
	open (OZ, $ZMOSF);
	<OZ>;
	while (<OZ>) {
		chomp;
		my ($osg, $ospep, @zmhittxts) = split("\t", $_);
		next if (@zmhittxts == 0);
		my $osclust = $refsp_pep_clust{"Oryza_sativa"}{$ospep}; 
		my @peps = @{ $refsp_clust_peps{"Oryza_sativa"}{$osclust}};
		next unless ($ospep eq $peps[0] ); # only get hits for top ranked (longest) pep in cluster 
		die "? file $ZMOSF has multiple entries for query ". $peps[0] . "\n" if (exists $osclustcheck{$osclust});
		$osclustcheck{$osclust} = 1;
		my ($zmhit, $hitdet) = split(" ", $zmhittxts[0]); # top scoring hit
		die "? zmclust not defined for $zmhit\n\n" unless (exists $refsp_pep_clust{"Zea_mays"}{$zmhit}); 	
		my $zmclust = $refsp_pep_clust{"Zea_mays"}{$zmhit}; 
		$hitdet =~ /score=(\d+)\;/; 
		my $score = $1;
		next if (exists $zmclust_tophit{$zmclust} && ($zmclust_tophit{$zmclust}{"score"} >= $score)); # ensure Zm hits are exclusive for each cluster so cannot map to multiple Os
		$zmclust_tophit{$zmclust} = {"score"=>$score, "osclust"=>$osclust}; 
	}
	close OZ;
	#+ reorganise to output structure
	my %osclust_blastp_zmclust;
	foreach my $zmclust (keys %zmclust_tophit) {
		$osclust_blastp_zmclust{$zmclust_tophit{$zmclust}{"osclust"} } = $zmclust;
	}
	return \%osclust_blastp_zmclust;
}


sub sort_orths {
	my %orth = %{ shift() };
	# return pep ids sorted by descending confidence, then percentage id to query seq
	my @ids = sort {$orth{$b}{confidence} <=> $orth{$a}{confidence}
					||
					$orth{$b}{q_pcid} <=> $orth{$a}{q_pcid}
					||
					$a cmp $b} keys %orth; # if all metrics equal, sort alphabetically for consistent output
	return \@ids;
}

sub set_orth_excsets {
# sets the topclust for each orth to restrict to this mapping and defines pepset of seeds that are excluded from ortholog set
	my $clustersetR = shift();
	#+ define topclust for each orth
	my %otc;
	foreach my $refsp (@REFSPS) {
		foreach my $clust (sort keys %{$refsp_clust_allorth{$refsp}} ) {
			next unless (exists $clustersetR->{$clust});
			foreach my $sp (sort keys %{$refsp_clust_allorth{$refsp}{$clust}} ) {
				foreach my $orth (sort keys %{$refsp_clust_allorth{$refsp}{$clust}{$sp}} ) {
					if (exists $otc{$orth} ) {
						my $topclust = $otc{$orth}{"clust"};
						my $tcrefsp = $otc{$orth}{"refsp"};
						$otc{$orth} = {"clust"=>$clust, "refsp"=>$refsp} if 
						  ($refsp_clust_allorth{$refsp}{$clust}{$sp}{$orth}{"t_pcid"} > $refsp_clust_allorth{$tcrefsp}{$topclust}{$sp}{$orth}{"t_pcid"});
					} else {
						$otc{$orth} = {"clust"=>$clust, "refsp"=>$refsp};
					}
				}
			}
		}
	}
	%orth_topclust = %otc; # re-set global variable
	#- 
	#+ define all cluster members from clusterset to be excluded as orths
	my %ps;
	foreach my $clust (keys %$clustersetR) {
		foreach my $refsp (@REFSPS) {
			if (exists $refsp_clust_peps{$refsp}{$clust}) {
				my @peps = @{ $refsp_clust_peps{$refsp}{$clust}};
				@ps{ @peps } = 1; # set hash slice
			}
		}
	}
	#-
	%pepset = %ps; # re-set global variable
	return;
}

sub get_exc_orth {
	# get orths which are exc i.e. do not match better to others in set
	# usage:
	#	my $excorthsR = Ortho_utils::get_exc_orth($refsp, $clust, $sp, $excseeds)
	my ($refsp, $clust, $sp, $excseeds) = @_;
	# following are accessed as global variables (large or unchanged):
	# %refsp_clust_allorth %orth_topclust
	#
	# output of sub:
	my %excorths;
	if (exists $refsp_clust_allorth{$refsp}{$clust}{$sp}) {
		foreach my $orth (sort keys %{ $refsp_clust_allorth{$refsp}{$clust}{$sp}} ) { 
			next if ($excseeds && exists $pepset{ $orth });
			my $topclust = $orth_topclust{$orth}{"clust"};
			# add orth if this is topclust or topclust is deleted from table so there is no overlap
			$excorths{$orth} = $refsp_clust_allorth{$refsp}{$clust}{$sp}{$orth} if ($clust eq $topclust ); 
		}
	}
	return \%excorths;
}

sub set_clustset_from_table {
	my %cs;
	foreach my $grpid (keys %grp_table) {
		my $refsp_clustR = _seedclust_from_grpid ($grpid);
		foreach my $c (values %$refsp_clustR) {
			$cs{ $c } = 1; 
		}
	}
	return \%cs;
}

sub add_orths_to_table {
	# takes grp table and adds exc orths for that table 
	# usage Ortho::add_orths_to_table(); #  used to exclude all cluster members 
	# following are accessed as global variables (large or unchanged)

	# get orths which are exc i.e. do not match better to others used in table and add ids to table
	foreach my $grpid (sort keys %grp_table) {
		my $refsp_clustR = _seedclust_from_grpid ($grpid);
		my @refsps;
		if (exists $refsp_clustR->{"Oryza_sativa"} && exists $refsp_clustR->{"Zea_mays"}) { # clust is seeded from both ref spp
			@refsps = ("Oryza_sativa", "Zea_mays");
		} elsif (exists $refsp_clustR->{"Oryza_sativa"}) { # clust is seeded from Os only
			@refsps = ("Oryza_sativa");
		} else {
			@refsps = ("Zea_mays");
		}
SP:		foreach my $sp (@sps) {
			my %refsp_orths; # exc orths for this grp and sp 
			foreach my $refsp (@refsps) {
				next SP if ($refsp eq $sp); # refsp are populated by cluster members (seeds) not orths
				my $clust = $refsp_clustR->{$refsp};
				$refsp_orths{$refsp} = get_exc_orth($refsp, $clust, $sp, 1);
				# foreach my $p (keys %{$refsp_orths{$refsp}} ) { # remove cluster members of other grps
					# delete $refsp_orths{$refsp}{$p} if (exists $pepset{ $p });
				# }
			}
			#	now add ortho ids to table	
			my @orthids;
			if (exists $refsp_clustR->{"Oryza_sativa"} && exists $refsp_clustR->{"Zea_mays"}) { # grp is seeded from both ref spp
				@orthids = @{ sort_orths( $refsp_orths{"Oryza_sativa"} )} if (exists $refsp_orths{"Oryza_sativa"});
				# get any extra orths from Zm blastp hit
				my @zmorthids = @{ sort_orths( $refsp_orths{"Zea_mays"} )} if (exists $refsp_orths{"Zea_mays"});
				foreach my $zmorthid (@zmorthids) {
					push @orthids, $zmorthid unless (exists $refsp_orths{"Oryza_sativa"}{$zmorthid});
				}
			} else { # grp is seeded from one ref sp only
				@orthids = @{ sort_orths( $refsp_orths{$refsps[0]} )} if (exists $refsp_orths{$refsps[0]});
			}
			
			$grp_table{$grpid}{$sp} = [@orthids] if (@orthids > 0);

		}# end each sp
	}

}

sub del_table_rows_below_min {
	# delete table entries with fewer than min num spp and return number deleted
	# usage:
	#	my $ndel = Ortho::del_table_rows_below_min($nspmin)
	my $nspmin = shift();
	my $ndel = 0;
	foreach my $grpid (sort keys %grp_table) {
		my $n = scalar keys %{$grp_table{$grpid}};
		if ($n < $nspmin) {
			delete $grp_table{$grpid};
			$ndel++;
		}
	}
	return $ndel;
}
	
sub _seedclust_from_grpid {
	my $grpid = shift();
	my %refsp_clust; 
	if ($grpid =~ /grp\:(\S*)(Zm\S*)$/) { # ** differs from v1.3 **
		$refsp_clust{"Oryza_sativa"} = $1 if (length $1>0);
		$refsp_clust{"Zea_mays"} = $2 if (length $2>0);
	} elsif ($grpid =~ /grp\:(\S*)$/) {
		$refsp_clust{"Oryza_sativa"} = $1 if (length $1>0);
	} else {
		die "cannot parse group_id $grpid\n";
	}
	die "no cluster ids found in $grpid\n" unless (exists $refsp_clust{"Oryza_sativa"} || exists $refsp_clust{"Zea_mays"});
	return \%refsp_clust;
}

sub nsp_stats {	# print summary stats for n clusts which have an entry for at least n spp
	# following are accessed as global variables (large or unchanged)
	# our @sps;
	# our %grp_table; 
	my %nsp = (9=>0, 12=>0, 15=>0);
	my $maxnsp = scalar @sps;
	$nsp{$maxnsp} = 0;
	foreach my $grpid (sort keys %grp_table) {
		my $n = scalar keys %{$grp_table{$grpid}};
		foreach my $nl (keys %nsp) {
			$nsp{$nl}++ if ($n >= $nl);
		}
	}
	foreach my $nl (sort {$a<=>$b} keys %nsp) {
		print $nsp{$nl} . " grps with a non blank entry >= $nl spp\n";
	}
	print "\n";
}

sub write_orth_table {
	# Usage:
	# my ($nout, $npotaltprof, $otf) = Ortho::write_orth_table();
	my $nout=0; # n groups output
	my $npotaltprof=0; #  number of potential alternative profiles to generate in next steps
	open (OTF, ">$ORTHOTABLEF");
	my @h = ("grp_id", @sps);
	print OTF join("\t", @h) . "\n";
	foreach my $grpn (sort keys %grp_table) {
		my @row = ($grpn);
		foreach my $sp (@sps) {
			if (exists $grp_table{$grpn}{$sp}) {
				my @peps = @{ $grp_table{$grpn}{$sp} };
				push @row, "@peps";
				$npotaltprof += scalar @peps - 1;
			} else {
				push @row, "";
			}
		}
		print OTF join("\t", @row) . "\n";
		$nout++;
	}
	close OTF;
	return ($nout, $npotaltprof, $ORTHOTABLEF);
}
	
sub read_orth_table {
	# Usage:
	# my ($OTF) = Ortho::read_orth_table();
#	my $expncol = scalar @sps + 1;
	open (OTF, $ORTHOTABLEF);
	my $hl = <OTF>;
	while (my $l = <OTF>) {
		chomp $l;
		my ($grpn, @cs) = split("\t", $l);
		my $grpid = $grpn;
		$grpid =~ s/grp\://; 
		for (my $i=0; $i<@cs; $i++) {
			if (defined $cs[$i] && length $cs[$i]>0) {
				my @peps = split(" ", $cs[$i]);
				$grp_table{$grpid}{$sps[$i]} = [@peps];
			}
		}
	}
	close OTF;
	return $ORTHOTABLEF;
}

1;
