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

my $hmmset = $ARGV[0]; # the name of HMM set just completed is passed so that a copy can be made

if (-e $FSSTATFALL) {
	print "$FSSTATFALL already exists.\n";
	my $cF = $FSSTATFALL . ".older";
	`mv -f $FSSTATFALL $cF`;
	if (-e $cF) {
		print "moved older $FSSTATFALL to $cF\n";
	} else {
		die "could not create $cF\n";
	}
}
		
my $h;
my %id_l;
foreach my $F (@FSSTATFS) {
	die "file $F not found\n" unless (-e $F);
	open (IN, $F);
	$h = <IN>;
	my $nl=0;
	while (my $l=<IN>) {
		my ($id, @a) = split("\t", $l);
		$id_l{$id} = $l;
		$nl++;
	}
	close IN;
	print $nl . " lines from $F\n";
}
my $nltot=0;
open (OUT, ">$FSSTATFALL");
print OUT $h;
foreach my $id (sort keys %id_l) {
	print OUT $id_l{$id};
	$nltot++;
}
close OUT;
print "Wrote $nltot lines to $FSSTATFALL\n";

if (defined $hmmset) {
	my $cpname = $FSSTATFALL;
	$cpname =~ s/$FSSTATFROOT/$hmmset/;
	`cp $FSSTATFALL $cpname`;
	print "made copy of status file to keep $hmmset info in:\n$cpname\n";
}
