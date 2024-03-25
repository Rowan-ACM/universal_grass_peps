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
# split status file for running in parallel


my %f_nl;
foreach my $F (@FSSTATFS) {
	$f_nl{$F} = 0;
	if (-e $F) {
		my $oF = $F . ".older";
		`cp -f $F $oF`;
		warn "file $F already exists. Made copy to $oF before overwrite\n"  ;
	}
}

open (IN, $FSSTATFALL);
my $h = <IN>;
my %id_l;
while (my $l=<IN>) {
	my ($id, @a) = split("\t", $l);
	$id_l{$id} = $l;
}
close IN;
my $nl = scalar keys %id_l;
print "$nl lines in $FSSTATFALL\n";

my $lastj = -1;
my $F;
my $nlperfile = 1 + $nl / @FSSTATFS;
print "approx $nlperfile lines per file\n";
my $i=0;
foreach my $l (values %id_l) { # this is just to randomise order of lines
	my $j = int ($i / $nlperfile); 
	if ($j > $lastj) {
		close OUT if (defined $F);
		$F =  $FSSTATFS[$j];
		open (OUT, ">$F") || die "cannot open $F\n";
		print OUT $h;
	}
	$lastj = $j;
	print OUT $l;
	$i++;
	$f_nl{$F}++;
}
close OUT;

foreach my $F (sort keys %f_nl) {
	print "wrote ". $f_nl{$F} . " lines to $F\n";
}