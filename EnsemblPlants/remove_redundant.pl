#!/usr/bin/env perl

##################################################
# universal_grass_peps pipeline v1.3             #
# Rowan Mitchell 2023                            #
# rowan.mitchell@rothamsted.ac.uk                #
# https://www.rowanmitchell-grassscience.co.uk/  #
##################################################

use strict;
use warnings;
use Bio::SeqIO;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;

my $SPF = "grass_sp_status.txt";
my $REDSEQF = "redundant_ids.txt"; # record of seqs that have been removed (used for lookup e.g. in ortho comparison)
my $GENESF = "gene_info.txt"; # gene info for all pep ids

my $gmfinished = Pipeline::read_gmsteps_status();
die "All gene model steps already completed." if ($gmfinished);

# removes identical redundant peptide seqs and outputs alphabetically lowest ID only of redundant set
open (OUT, ">$GMDIR$REDSEQF");
open (GEN, ">$GMDIR$GENESF");
foreach my $sp (@sps) {
	print "**$sp** doing..\n";
	my $fF =$sp_stat{$sp}{ori_pep_fasta_file};
	my $outF =$fF;
	$outF =~ s/\.all\./\.nr\./;
	my $n=0;
	
	my $in = new Bio::SeqIO( -format => 'fasta', -file   => "$GMDIR$fF" );
	my %seq_ids;
	my %id_seqobj;
	while (my $seqobj = $in->next_seq() ) {
		my $id = $seqobj->display_id();
		$id_seqobj{ $id } = $seqobj;
		my $seq = $seqobj->seq();
		push @{ $seq_ids{ $seq}}, $id;
		$n++;
		my $desc = $seqobj->desc();
		my $gene;
		if ($desc =~ / gene\:(\S+) /) {
			$gene = $1;
		} else {
			$gene = "gene not found in $desc";
		}
		print GEN "$sp\t$id\t$gene\n";
	}
	print "read in $n seqs from $fF\n";

	my $nout=0;
	my $out = new Bio::SeqIO( -format => 'fasta', -file   => ">$GMDIR$outF" );
	foreach my $seq (keys %seq_ids) {
		my @ids = sort @{ $seq_ids{ $seq}};
		if (@ids > 1) {
			for (my $i=1; $i<@ids; $i++) {
				print OUT "$sp\t$ids[$i]\t$ids[0]\n";
			}
		}
		my $seqobj = $id_seqobj{ $ids[0] };
		$out->write_seq( $seqobj ); 
		$nout++;
	}

	print "\tInput $n seqs. Output $nout unique peptide seqs to $outF\n";

}
close GEN;
close OUT;
