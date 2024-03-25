#!/usr/bin/env perl

##################################################
# universal_grass_peps pipeline v1.4             #
# Rowan Mitchell 2020-2024                       #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################


use strict;
use warnings;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;
use Ortho;

# get spp
my ($gmfinished, $spsR) = Pipeline::read_gmsteps_status();	
#
# following are passed to subroutines as global variables (large or unchanged). Exported from Orth
# our %refsp_clust_peps; # clustering of peps into groups from self blastp clust -> [@peps]
# our %refsp_pep_clust; # clustering of peps into groups from self blastp	pep -> clust
# our %refsp_clust_allorth; # all orthologs for each refsp clust
# our %grp_table; # Table of orths that is main output of script

my %clusterset; # set of clusters to include for defining orth_topclust
foreach my $refsp (@REFSPS) {
	print "\nRef species $refsp\n";
	
	my $nc;
	
	($refsp_clust_peps{$refsp}, $refsp_pep_clust{$refsp}) = Ortho::get_clusters($refsp); # all member peps from  grps as hash clust -> [@members]
	$nc = keys %{ $refsp_clust_peps{$refsp}};
	print "read in $nc clusters from self blastp for $refsp\n";

	#+ inputs ortholog tables to find potential sets of genes in all grasses
	my $orthoF = $REFSP_SHTN{$refsp} . $OUTORTHOEXT;
	open (IN, $orthoF) || die "cannot open $orthoF\n";
	<IN>;
	while (<IN>) {
		chomp;
		my ($clust, @cells) = split("\t", $_);
		for (my $i=0; $i<@cells; $i++) {
			my $sp = $sps[$i];
			my @cs = split(";", $cells[$i]);
			foreach my $c (@cs) {
				my $orth;
				my %orthdet;
				($orth, @orthdet{@ORTHKEYS}) = split(" ", $c);
				$refsp_clust_allorth{$refsp}{$clust}{$sp}{$orth} = \%orthdet;
				$clusterset{$clust}++;
			}
		}
	}
	close IN;
	#- inputs ortholog tables to find potential sets of genes in all grasses
	$nc = keys %{ $refsp_clust_allorth{$refsp}};
	print "Read in all orths from $orthoF . $nc $refsp clusters have at least 1 orth\n";
	
}

#+ get Os -> Zm blastp mapping
my $osclust_blastp_zmclustR = Ortho::read_in_Zm_blastp_Os();
print "\nGot Zm - Os blastp mapping from $ZMOSF\n";
#-

#+ --------------------- start to define maximal set of groups -----------
# %grp_table (imported from Ortho) # hash keyed by grp name -> table of spp. Each sp -> array of pep ids in descending priority 
	#   grp name ($grpn) is "grp:" first Os pep (if any) . first Zm pep (if any). So one of these refsp must be present
	#  below defines seed grps and adds all peps for these from clusters 
	# When these are defined. orthologs are added in subroutine

my %zmclust_notos;
@zmclust_notos {keys %{ $refsp_clust_allorth{"Zea_mays"} }} = 1; # nr list of all Zm clusters with >0 orths
# this is used for rare cases where no good Os gene model so needs to be found by genblast in Os genome

Ortho::set_orth_excsets(\%clusterset); #define top cluster for each orth

my $nocz=0; my $noczbp=0;
foreach my $osclust (sort keys %{$refsp_clust_allorth{"Oryza_sativa"}} ) {

	my $oszmorthsR = Ortho::get_exc_orth("Oryza_sativa", $osclust, "Zea_mays", 0); # find any exc Zm orths for this Os clust

	my $grpn;
	if (keys %$oszmorthsR > 0) { # has Zm orths 
		$grpn = "grp:$osclust";
		
		##+ following just for stats output
		$nocz++;
		my @zmorthids = @{ Ortho::sort_orths( $oszmorthsR ) }; 
		# delete all zmclust_notos clusters which have orthologs to Os
		foreach my $zmorthid (@zmorthids) {
			my $zmclust = $refsp_pep_clust{"Zea_mays"}{$zmorthid};
			delete $zmclust_notos{$zmclust} if (exists $zmclust_notos{$zmclust});
		}
		##- 
		
	} elsif (exists $osclust_blastp_zmclustR->{$osclust}) {			# has no Zm orths but does have Zm blastp hits
		my $zmclust = $osclust_blastp_zmclustR->{$osclust}; 	
		die "Zm cluster id $zmclust does not start with \'Zm\'. Needed for parsing\n" unless ($zmclust =~ /^Zm/);
		$grpn = "grp:$osclust$zmclust"; # this indicates table grp is seeded by both Os and Zm to find orthologs ** differs from v1.3 **
		$grp_table{$grpn}{"Zea_mays"} = $refsp_clust_peps{"Zea_mays"}{$zmclust} ; # copy array ref for seed Zm peps

		##+ following just for stats output 
		$noczbp++;
		delete $zmclust_notos{$zmclust} if (exists $zmclust_notos{$zmclust});# delete zmclust_notos self grps found by blastp from Os
		##- 

	} else {
		$grpn = "grp:$osclust"; #  "No Zm orths or blastp hits for Os grp $osclust\n";
	}
	
	$grp_table{$grpn}{"Oryza_sativa"} = $refsp_clust_peps{"Oryza_sativa"}{$osclust} ;  # copy array ref for seed peps

}
my $noc = keys %{$refsp_clust_allorth{"Oryza_sativa"}};
print "\nOut of $noc Oryza_sativa clusters $nocz have Zm orthos and $noczbp have no Zm orthos but have Zm blastp hits\n";
my $nzc = keys %{ $refsp_clust_allorth{"Zea_mays"} };
my $nzc_notos = keys %zmclust_notos;
print "\nOut of $nzc Zea_mays clusters $nzc_notos are not orthologous or similar to Os. Adding grps for these\n";

# add entries for Zm seeds only
foreach my $zmclust (sort keys %zmclust_notos ) {

	my $grpn = "grp:$zmclust"; # Zm ids are recognised by pattern "Zm"
	$grp_table{$grpn}{"Zea_mays"} = $refsp_clust_peps{"Zea_mays"}{$zmclust} ; # copy array ref for seed peps

}
#- max list grps now defined
my $nctot = keys %clusterset;
%clusterset = %{ Ortho::set_clustset_from_table() }; # redefine clusterset to those in table
Ortho::set_orth_excsets(\%clusterset); #define top cluster for each orth

my $nctab = keys %clusterset;
print "\nNumber clusters: All $nctot. Now in table $nctab.\n";

Ortho::add_orths_to_table();
print "\nTable with Os and Zm orths plus Zm only orths:\n";
Ortho::nsp_stats();

#+ now delete table entries below min num spp and reassign exc orthos
my $nspmin = scalar @sps - $NG_MAX_MISSGM;
my $ndel = Ortho::del_table_rows_below_min($nspmin);
%clusterset = %{ Ortho::set_clustset_from_table() }; # redefine clusterset to those in revised table
Ortho::set_orth_excsets(\%clusterset); #define top cluster for each orth
$nctab = keys %clusterset;
print "$ndel entries in table removed as below min num spp $nspmin. Number clusters now in table $nctab.\n";
Ortho::add_orths_to_table(); # re-do orths to allow ones from del rows to be reassigned (only small number affected)
print "\nTable with Os and Zm orths plus Zm only orths. Entries below $nspmin spp removed:\n";
Ortho::nsp_stats();
#-

#main output
my ($nout, $npotaltprof, $otf) = Ortho::write_orth_table();
print "Wrote $nout grps with at least $nspmin spp to $otf\n";
print "profile refinement for HMM-0 could require generation of $npotaltprof profiles (maximum)\n";

exit(0);
