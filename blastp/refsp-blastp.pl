#!/usr/bin/env perl

##################################################
# universal_grass_peps pipeline v1.3             #
# Rowan Mitchell 2023                            #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################


# 1. generates blastp results for reference sp. : Os - Os, Zm -Zm and Zm - Os
# 2. from self blastp clusters peps that are >90% id. Cluster name is just longest pep id
#    this is output together with all single peps to files $refsp."-self_groups.txt"
# 3. output all exc hits for every Os from Zm queries to "Zea_mays-Os.hits.txt"


use strict;
use warnings;
use Bio::SeqIO;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;
use Blastp;

#+ these are kept from loop for Os as needed for Zm -> Os
my %osid_len;
my %ospep_osclustern; # #key is all nr Os pep, value is os clustername
my %osclustern_ospeps; # reverse hash 
#-

foreach my $refsp (@REFSPS) {
	print "\nDoing self blastp for $refsp\n\n";
	my $faF = Pipeline::gmfasta_filename( $refsp ); # gene model pep fasta filename
	my $blastdbF = $faF . ".pot";
	if (-e $blastdbF) {
		print "blastdb for $refsp already exists.\n";
	} else {
		print "making blastdb for $refsp\n";
		`makeblastdb -in $faF -dbtype prot`;
	}
	print "running self blastp for $refsp\n";
	my ($blastpF, $nql) = Blastp::run_blastp($refsp, $refsp); # self blast
	print "$nql queries in blastp result file $blastpF\n";

	#------------------------------- analyse blastp file -----------------------------
	
	print "\nAnalysing self blast for ref sp $refsp\n";
	my ($ids_ref, $id_len_ref) = Pipeline::read_fasta_lens( $faF );
	
	my $hits_ref = Blastp::get_blastp_hits($blastpF, $id_len_ref);
	my %hits = %{ $hits_ref};
	#-----------

	my %memb_cluster; #key is member pep id, value is cluster_id
	my %cluster_membs; # key is cluster_id, value is %members
	my $npep=0;
	foreach my $q (@$ids_ref) {
		unless (exists $hits{$q}{$q}) {
			print "$q\tno self hit?\n";
			next;
		}
		$npep++;
		foreach my $s (keys %{ $hits{$q}} ) {
			next if ($q eq $s);
			next unless (exists $hits{$s}{$s});
	# to calculate relative score use bigger ( usually longer) out of query and subject self-score
	#		my $s_ss = $hits{$s}{$s}{score}; # subject self score
	#		my $fs;

			next unless (exists $hits{$s}{$q}); # reciprocal search should find query if good hit
			my $pwpc_id; # pairwise %id set to lower of reciprocal searches
			if ($hits{$q}{$s}{pc_id} > $hits{$s}{$q}{pc_id}) {
				$pwpc_id = $hits{$s}{$q}{pc_id};
			} else {
				$pwpc_id = $hits{$q}{$s}{pc_id};
			}
			next if ($pwpc_id < $SELFBLASTPLIM);
			my %h = ($q=>$pwpc_id, $s=>$pwpc_id);
			if (exists $memb_cluster { $q}) {
				die "corresponding cluster_membs doesn't exist for key $memb_cluster{$q}\n" 
				if (! exists $cluster_membs{ $memb_cluster { $q} });
				my %q_cluster = %{ $cluster_membs{ $memb_cluster { $q} }}; 
				foreach my $q_mem (keys %q_cluster) {	# add all other members of q group to %h
					next if (exists $h{ $q_mem } && $h{ $q_mem } > $q_cluster{ $q_mem });
					$h{ $q_mem } = $q_cluster{ $q_mem };
				}
			}
			if (exists $memb_cluster { $s}) {
				die "corresponding cluster_membs doesn't exist for key $memb_cluster{$s}\n" 
				if (! exists $cluster_membs{ $memb_cluster { $s} });
				my %s_cluster = %{ $cluster_membs{ $memb_cluster { $s} }}; 
				foreach my $s_mem (keys %s_cluster) {   # add all other members of s group to %h
					next if (exists $h{ $s_mem } && $h{ $s_mem } > $s_cluster{ $s_mem });
					$h{ $s_mem } = $s_cluster{ $s_mem };
				}
			}
			#+ find/set group id for members in %h
			my $cluster_id; # numerical id for cluster 
			if (exists $memb_cluster { $q} && exists $memb_cluster { $s}) { # both belong to existing clusters
				if ($memb_cluster { $q} > $memb_cluster { $s}) { # if different, merge and set to lower id
					delete $cluster_membs{ $memb_cluster { $q} }; 
					$cluster_id = $memb_cluster { $s};
				} elsif ($memb_cluster { $q} < $memb_cluster { $s}) { # if different, merge and set to lower id
					delete $cluster_membs{ $memb_cluster { $s} }; 
					$cluster_id = $memb_cluster { $q};
				} else { # same group
					$cluster_id = $memb_cluster { $q};
				}
			} elsif (exists $memb_cluster { $q}) { # only 1 exists, use this id
				$cluster_id = $memb_cluster { $q};
			} elsif (exists $memb_cluster { $s}) { # only 1 exists, use this id
				$cluster_id = $memb_cluster { $s};
			} else { # neither in cluster, create new one
				my @cluster_ids = sort {$a <=> $b} keys %cluster_membs;
				if (@cluster_ids eq 0) {
					$cluster_id = 1;
				} else {
					$cluster_id = $cluster_ids[$#cluster_ids] + 1; # new cluster id is 1 more than last element of sorted keys 
				}
			}
			#-
			# add %h to clusters and set members in reciprocal hash to correct group id
			$cluster_membs{ $cluster_id }  = \%h;
			foreach my $memb (keys %h) {
				$memb_cluster { $memb} = $cluster_id;
			}
		}
	}
		

	# get stats
	my $size_biggest=0;
	my $ncluster = keys %cluster_membs;
	my $nmemb=0;

	my %pep_clustern; # put %cluster_membs in new structure for all nr peps : first member -> array all members
	my %clustern_peps; # reverse hash 
	foreach my $cluster (keys %cluster_membs) {
		# sort member peps by descending length (then by id alphabetically)
		my @members = sort 
			{$id_len_ref->{$b} <=> $id_len_ref->{$a} || $a cmp $b} 
			keys %{ $cluster_membs {$cluster} };
		my $clustern = "cluster:" . $members[0]; # group name is made from first (longest) pep
		foreach my $pep (@members) {$pep_clustern{ $pep} = $clustern;}		
		$clustern_peps{$clustern} = [ @members ]; 
		$size_biggest = @members if (@members > $size_biggest);
		$nmemb += @members;
	}
	# self clusters stats output
	my $pc_cluster = $nmemb / $npep * 100.;
	print "\nWith limit on \% id of $SELFBLASTPLIM\% :\n";
	print "$nmemb out of $npep pep are in $ncluster groups\n";
	print "i.e. $pc_cluster \% in groups\n";
	print "biggest cluster has $size_biggest members\n";
	
	# now add single member clusters
	foreach my $id (@$ids_ref) {
		unless (exists $pep_clustern{$id}) {
			my $clustern = "cluster:" . $id; # group name is made from pep id
			$pep_clustern{ $id} = $clustern;
			$clustern_peps{ $clustern } = [$id];
		}
	}
	
	if ($refsp eq "Oryza_sativa") {# keep data for Zm->Os later
		%osid_len =%{ $id_len_ref} ; 
		%ospep_osclustern = %pep_clustern; 
		%osclustern_ospeps = %clustern_peps; 
	}

	# OUTPUT
	my $outF = $refsp . $CLUSTFEXT;
	open (OUT, ">$outF");
	my @h = ("cluster_id");
	for (my $i=1; $i<=$size_biggest; $i++) {
		push @h, "pep$i";
	}
	print OUT join("\t", @h) . "\n";
	foreach my $clustern (sort keys %clustern_peps) {
		print OUT join("\t", ($clustern, @{ $clustern_peps{$clustern} })) . "\n";
	}
	close OUT;
	print "Wrote $outF\n";
}

print "\nDoing Zm->Os steps\n\n";
my ($zmosblastpF, $zmosnql) = Blastp::run_blastp("Zea_mays", "Oryza_sativa"); # Zm -> Os blast
print "$zmosnql queries in blastp result file $zmosblastpF\n";
my $zmoshits_ref = Blastp::get_blastp_hits($zmosblastpF, \%osid_len);
my %oscluster_exchits = %{  Blastp::get_subject_sets_exchits ($zmoshits_ref, \%ospep_osclustern) };

open (OUT, ">$ZMOSF"); # main output file - exclusive hits
print OUT join("\t", "group_id", "Os pep", "Zea_mays hit1", "Zea_mays hit2...") . "\n";
foreach my $clustern (sort keys %oscluster_exchits) {	
	foreach my $ospep ( @{ $osclustern_ospeps{ $clustern }} ) {
		my @row = ( $clustern, $ospep );
		if (exists $oscluster_exchits{$clustern}{$ospep}) {
			foreach my $zmpep (@ {$oscluster_exchits{$clustern}{$ospep}}) {
				my $det = Blastp::write_hitdetail( $zmoshits_ref->{$zmpep}{$ospep} );			
				push @row, "$zmpep $det";
			}
		}
		print OUT join("\t", @row) . "\n";
	}
}
close OUT;
print "Output written to $ZMOSF \n";

exit(0);
