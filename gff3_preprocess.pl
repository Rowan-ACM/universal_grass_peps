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
use Bio::Tools::GFF;
use lib "/home/data/tropical_forage_grasses/grass_core/final/perl_lib";
use Pipeline;
use Genblast;

my %GFFPRIMARYTAGS = ("ensembl"=>{"cds"=>"mRNA", "exon"=>"CDS"},
                      "genblast"=>{"cds"=>"transcript", "exon"=>"coding_exon"});

my $sp = $ARGV[0] || "Brachypodium_distachyon";

my $gmfinished = Pipeline::read_gmsteps_status();	
die "gff3 filename not found for $sp\n" unless (exists $sp_stat{$sp}{gff3_file}); 
my $gffF = $GFFDIR . $sp_stat{$sp}{gff3_file};
die "gff3 file $gffF for $sp not there\n\n" unless (-e $gffF); 

my $gffppF = $gffF . "_pp.txt";


print "reading Ensembl file\n";	
my ($gsR, $geneR) = Genblast::parse_gff3("ensembl", $gffF);
my %gene = %{$geneR};

my $ng = keys %gene; 
print "$ng genes in $gffF\n";

my $chr;
my @gffppFs;
foreach my $g (@$gsR) {
	unless (exists $gene{$g}{exons} && exists $gene{$g}{cds}) {
		warn "gene $g lacks exons or cds\n";
		next;
	}
	my $gdtxt = Genblast::gene_details_txt($gene{$g});
	my @ds = split(";", $gdtxt);
	if (! defined $chr || $chr ne $ds[0]) {
		close OUT if (defined $chr);
		$chr = $ds[0];
		my $gffppF = $gffF . "\.$chr\_pp.txt";
		push @gffppFs, $gffppF;
		open (OUT, ">$gffppF");
	}	
	my $id = $g; $id =~ s\transcript:\\i;
	print OUT join("\t", ($gdtxt, $id)) . "\n";
}
close OUT;
print "Output written to @gffppFs \n";

exit(0);




