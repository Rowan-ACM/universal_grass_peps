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
use HMMs;
use Genblast;

# completes group steps for grassgenes pipeline when best groups have been estimated for existing Ensembl gene models in current annotations. 
# Role is to try genblastn for all missing or low scoring members to try to get as many groups as possible in final set that pass min score criterion

#+ steps to complete for each grp
# 1 find set to try genblast 
# 2 try genblast
# 3 add/substitute genblast peps to group 
# 4 re-create hmm profile
# 5 if no improvement, revert to original gene model
# 6 create final hmm 
# 7 when all done, output some stats
#-
#+ outputs 
# 1. final grp hmms and info on their composition
# 2. info on genblast improvements to gene models
#-

my $statfileopt = $ARGV[0] || die "supply which status file as first argument (1-12) or 'all' or subset filename\n";

# per sp status file (GM steps ; must be all complete). defines @sps
my $gmfinished = Pipeline::read_gmsteps_status();	
die "All genemodel need to be complete to run this\n" unless ($gmfinished);

# get grpid status
# %grpid_status holds info on how far along pipeline each grpid has progressed
my	($grpid_statusR, $fsstatF);  
my $usesubset = 0;
my $subsetF;
my $doall = 0;
if ($statfileopt =~ /^\d{1,2}$/ && $statfileopt >= 1 && $statfileopt <= 12) {
	($grpid_statusR, $fsstatF) = Pipeline::get_fsteps_status($statfileopt);
} else {
	($grpid_statusR, $fsstatF) = Pipeline::get_fsteps_status("all"); # all for all and subset
	$doall = ($statfileopt eq "all");
	$usesubset = (! $doall);
	$subsetF = $statfileopt if ($usesubset); # e.g. subset "genes_with_evidence.txt";
}

# find how many to do 
my @setgrpids; # dependent on subset or all
if ($usesubset) {
	my ($ssgrpidsR, $ssdb) = Pipeline::read_subset($subsetF, $grpid_statusR); 
	@setgrpids = @$ssgrpidsR;
} else { # all or set1-12
	@setgrpids = keys %$grpid_statusR;
}
my %grpidstodo;
@grpidstodo{@setgrpids} = 1;
foreach my $grpid (@setgrpids) {
	die "Step incorrect for $grpid :" . $grpid_statusR->{$grpid}{"step_complete"} . " only \'HMM2\' or \'genblast\' allowed\n" unless
    	($grpid_statusR->{$grpid}{"step_complete"} eq "HMM2" || $grpid_statusR->{$grpid}{"step_complete"} eq "genblast" );
	delete $grpidstodo{$grpid} if ($grpid_statusR->{$grpid}{"step_complete"} eq "genblast"); # already done
}
my $ng = keys %grpidstodo;
print "Using status file $fsstatF\nThere are $ng grpids to do to create final HMM trying genblast from $fsstatF\n";
 
# following are for some final statistics
my %genblast_ngm = (failed=>0, not_improved=>0, improved=>0); 
my %summ_ngrpid = (failed=>0, passed=>0); 
my @passgrpids;

# main loop
foreach my $grpid (sort keys %grpidstodo) {
	my $grpreporttxt = "\n********************* $grpid *******************************\n";
	
	die "expected group details not found in $fsstatF for $grpid\n" unless (exists $grpid_statusR->{$grpid}{"current_hmm_details"});
	my $grpdetailsR = $grpid_statusR->{$grpid}{"current_hmm_details"};
	
	my @missing;
	my @outliers;
	#+ --- check for genblast already done for sp as indicated by id name	
	my $hmmFsR = Pipeline::group_fnames($grpid);
	my $id_spR = Pipeline::get_sp_from_faF($hmmFsR->{fasta});
	my $hasgenblast=0;
	my @gbids;
	foreach my $id (keys %$id_spR) {
		if ($id =~ /genblast/) {
			$hasgenblast = 1;
			push @gbids, $id;
		}
	}
	warn "? genblast ids already used in HMM for $grpid: @gbids\n" if ($hasgenblast);	
	#-
		
	my %members = %{ $grpdetailsR->{members} };
	# find missing and outliers
	foreach my $sp (@sps) { 
		push @missing, $sp if (! exists $members{$sp}); 
	}
	my $ntoadd = $NG_MAX_MISSGM - scalar @missing; # add up to this number more seqs to retry if score too low
	my @sorted_sps = sort {$members{$a}{hrss} <=> $members{$b}{hrss}} keys %members;
	my $devtxt = "";
	for (my $i=0; $i<$ntoadd; $i++) {
		my $sp = $sorted_sps[$i];
		if ($members{$sp}{hrss} < $MINRELSC || ($members{$sp}{hrss} < 0.9 && $members{$sp}{hrss_dev} < $DEVLIM) ) {
			push @outliers, $sp;
			$devtxt .= sprintf "%.3f %.1f ", $members{$sp}{hrss}, $members{$sp}{hrss_dev};
		}
	}
	$grpreporttxt .= "\t $grpid - missing @missing outliers @outliers $devtxt\n";
	# combine missing and outliers into %spstoredo with hrss values
	my %spstoredo;
	foreach my $sp (@missing) {
		$spstoredo{$sp} = 0;
	}
	foreach my $sp (@outliers) {
		$spstoredo{$sp} = $members{$sp}{hrss}; # set to initial score to check if improved
	}

	#-------- using genblast where needed
	my $ntoredo = scalar keys %spstoredo;    
	if ($ntoredo == 0) {
		$grpreporttxt .= "\tgenblast step. No missing or outliers so not needed\n";
	} else {
		$grpreporttxt .= "\tgenblast step. Trying genblast to improve $ntoredo gene models\n";
		my @genblastinfo; # info for status output
		my $changed = 0; # flag for at least 1 successful genblast 
		my $skiprest = 0; # flag to give up on genblast 
		my @sorted_spstoredo = sort {$spstoredo{$a} <=> $spstoredo{$b}} keys %spstoredo; # do missing first as if these fail abandon group
		foreach my $sp (@sorted_spstoredo) {
			my ($failed, $gbtxt, $gbpepfaF, $gbpepid) = Genblast::run_genblast( $sp, $grpid);
			if ($failed) { # failed means genblast found no gene model at all 
				$genblast_ngm{failed}++; 
				push @genblastinfo, sprintf "%s %.2f %s", "$sp ", $spstoredo{$sp}, "NO GENBLAST HIT";
				$grpreporttxt .= "\tGenblast: NO gene model found in $sp\n$gbtxt";
			} else {
				$grpreporttxt .= "\tGenblast for species $sp:\n$gbtxt\tPep outputfile created $gbpepid\n";
				my ($isbetter, $scorecmptxt, $newgrpdetailsR) = HMMs::try_replacing_seq_in_group($grpid, $grpdetailsR, $sp, $gbpepfaF, $gbpepid);
				if ($isbetter) {
					$grpdetailsR = $newgrpdetailsR;
					%members = %{ $grpdetailsR->{members} };
					$grpreporttxt .= "\tGenblast: success! improved group with $gbpepid\: $scorecmptxt\n";
					$changed++;
					$genblast_ngm{improved}++; 
					Genblast::keep_genblast_files($grpid);
				} else {
					$grpreporttxt .= "\tGenblast: NO improvement group with $gbpepid\: $scorecmptxt\n";
					$genblast_ngm{not_improved}++;  
				}
				push @genblastinfo, $scorecmptxt;
			}
			Genblast::remove_genblast_files($grpid); # remove all output on genblast dir for this grp
			my $currentscore = 0;
			$currentscore = $members{$sp}{hrss} if (exists $members{$sp});
			$skiprest = ($currentscore < $MINRELSC); # if score is still below min for this sp (will be 0 if still missing), no point in doing others
			if ($skiprest) {
				$grpreporttxt .= "\t$grpid score $sp " . $currentscore . " is too low. Skipping any remaining genblast\n";
				push @genblastinfo, sprintf "%s%.2f%s", " score below min (", $MINRELSC, ")";
				last;
			}
		}
		# 
#		#- 
		$grpreporttxt .= "\t$ntoredo tried with genblast. sp old_score new_score:". join(":", @genblastinfo) ."\n";
	}
	$grpid_statusR->{ $grpid}{"step_complete"} = "genblast";
	$grpid_statusR->{ $grpid}{"current_hmm"} = "final_hmm";
	$grpid_statusR->{ $grpid}{"current_hmm_details"} = $grpdetailsR;
	Pipeline::update_fsteps_status($grpid_statusR, $fsstatF);
# #--------final hmm file versions are present. get status text and check if pass
	my ($pass, $stattxt) = HMMs::final_hmmstatus($MINRELSC, $grpdetailsR);
	$grpreporttxt .= "\tgroup complete $grpid. Final status: $stattxt\n";
	if ($pass) {
		push (@passgrpids, $grpid);
		$summ_ngrpid{passed}++;
	} else {
		$summ_ngrpid{failed}++;
	}
#--------
	print $grpreporttxt; 
}

my $endreporttxt = "\n---------------------------------------------------------------------------\n";
$endreporttxt .= "\nGroup steps complete for set $statfileopt groups.\n";
$endreporttxt .= "See $fsstatF\n\n";
$endreporttxt .= "\n---------------------------------------------------------------------------\n";

$endreporttxt .= "\n********************* summary *******************************\n";
$endreporttxt .= "Genblast results - numbers of gene models:\n";
foreach my $k (sort keys %genblast_ngm) {
	$endreporttxt .= "\t $k ". $genblast_ngm{$k} . "\n";
}
$endreporttxt .= "group results with min score $MINRELSC- numbers of groups:\n";
foreach my $k (sort keys %summ_ngrpid) {
	$endreporttxt .= "\t $k ". $summ_ngrpid{$k} . "\n";
}
print $endreporttxt;	
exit(0);

