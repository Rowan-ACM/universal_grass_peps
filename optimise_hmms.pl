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
use Ortho;
use HMMs;

## Optimises grp membership to maximise hmm score at 2 different steps in pipeline
## first is at start to create first hmm HMM0 then refine this with alternative orthologs to HMM1
## second is to try all other non-member gene models. After this step, hmms are HMM2

my %STEPTODO_STEPDONE = ("HMM1"=>"start", "HMM2"=>"HMM1");

#+ start -> HMM1
# grass groups created generate-ortho-table.pl used as starting point
# for each of these groups, create a profile with first pep for each sp (HMM0)
# Then refine this by looking at scores of alt peps for each sp. against HMM0 If score is better than member by min amount, try substituting
# If member score and lowest score better, accept substitution.
# When refinement finished, resulting profile is HMM1
#-

#+ HMM1 -> HMM2
# grass groups created by this script in HMM1 mode used as starting point
# Then optimise this by looking at scores of alt peps (any non-member) for each sp. against HMM1 If score is better than member by min amount, try substituting
# If member score and lowest score better, accept substitution.
# When refinement finished, resulting profile is HMM2
#-

my $steptodo = $ARGV[0] || "none"; # step option from cl
my @ks = sort keys %STEPTODO_STEPDONE;
die "first argument $steptodo not recognised. Must be one of @ks\n" unless (exists $STEPTODO_STEPDONE{$steptodo});
my $statfileopt = $ARGV[1] || "all"; # status file option from cl

# per sp status file (GM steps ; should be all complete). defines @sps
my $gmfinished = Pipeline::read_gmsteps_status();	
die "All genemodel need to be complete to run this\n" unless ($gmfinished);

# get grpid status
my	($grpid_statusR, $fsstatF) = Pipeline::get_fsteps_status($statfileopt); 
print "Using status file $fsstatF\n";
my %grpidstodo;
foreach my $grpid (keys %$grpid_statusR) {
	if ($grpid_statusR->{ $grpid}{"step_complete"} eq $STEPTODO_STEPDONE{$steptodo}) { # stat is at previous step 
		$grpidstodo{$grpid} = 1; # so add to list to do
	} elsif ($grpid_statusR->{ $grpid}{"step_complete"} ne $steptodo) { # poss problem as stat is not equal to only other allowed value (current step done)
		warn "stat value of grp $grpid \'" . $grpid_statusR->{ $grpid}{"step_complete"} . "\' is unusual for step $steptodo. Expected to be $steptodo or "
		. $STEPTODO_STEPDONE{$steptodo} ."\n";
	}
}
my $ng = keys %grpidstodo;
print "There are $ng grpids to do to generate $steptodo from $fsstatF\n";

my $grpid_sp_difscoresR; # ONLY FOR HMM2 MODE
if ($steptodo eq "HMM1") {
	# get all ortholog table entries
	my ($OTF) = Ortho::read_orth_table();
	my $north = keys %grp_table;
	print "There are $north grpids defined by orthos in $OTF\n";
} else {	#HMM2
	$grpid_sp_difscoresR = HMMs::input_nonmem_scores($HMM1HITSF);# input scores for alt peps already calculated so just read in
}


my $nb = 0;
foreach my $grpid (sort keys %grpidstodo) {

	my $faFsR = Pipeline::group_fasta_fnames($grpid); # get fasta filenames 
	my $hmmFsR = Pipeline::group_fnames($grpid); # get  HMM filenames 

	#+ variables obtained differently according to whether HMM1 or HMM2 mode
	my $grpdetailsR; # current hmm details to be optimised
	my $altfaF; # file containing alt peps fasta
	my %altid_sp;	#  sp of alt peps 
	my %altid_scoredif; # holds alt ids which score better than original members
	#-
	
	if ($steptodo eq "HMM1") {
	#+++++++++++++++++++++ HMM0 PART. ONLY FOR HMM1 MODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		#+ check no existing files if doing from start
		my $hmmdone = 0;
		foreach my $k (sort keys %$hmmFsR) {
			if (-e $hmmFsR->{$k}) {
				my $hmmF = $hmmFsR->{$k};
				$hmmdone++;
			}
		}
		warn "$hmmdone hmm files already exist for $grpid . Will be overwritten\n" if ($hmmdone);
		#-

		my $hmm0faF = $faFsR->{"HMM0"};	
		#+ check HMM0 fa content is exactly as expected and assign sp for each alt id
		my %id_sp = %{ Pipeline::get_sp_from_faF($hmm0faF) };
		my %sp_id;
		while (my ($k, $v) = each %id_sp) {$sp_id{$v} = $k;}
		foreach my $sp (@sps) {
			if (exists $grp_table{$grpid}{$sp}) {
				my ($id, @altids) = @{ $grp_table{$grpid}{$sp} };
				die "expected $sp in grp $grpid not found\n" if (! exists $sp_id{$sp});
				die "expected $sp is in grp $grpid but id is $sp_id{$sp}\n" unless ($sp_id{$sp} eq $id);
				@altid_sp{@altids} = ($sp) x @altids;
			}
		}
		
		print "$hmm0faF composition ok\n";
		#- 
		#+ --------------  run hmmscan for alt orthos versus HMM0 db
		#--------------  create HMM0 database and get new selfscores
		$grpdetailsR = HMMs::make_groupdb($grpid, $hmm0faF);
		# [HMM0 step is quick and status not currently stored- always redone after interuption]
		print "HMM0_complete $grpid\n";
		$altfaF = $faFsR->{"HMM0alt"};	
		my $altid_scoreR = HMMs::get_alt_scores_groupdb($grpid, $altfaF); # # calc scores for alt peps 
		#+ ------ find alt ids which score better than original members
		foreach my $altid (keys %$altid_scoreR) {
			die "no sp found for alt id $altid\n" unless (exists $altid_sp{$altid} && defined $altid_sp{$altid} );
			my $sp = $altid_sp{$altid};
			my $scoredif = $altid_scoreR->{$altid} / $grpdetailsR->{maxposs_score} - $grpdetailsR->{"members"}{$sp}{"hrss"};
			$altid_scoredif{$altid} = $scoredif if ($scoredif > $MINHRSSDIF);
		} 
		#- ------ 

	#-----------------------------------------------------------------------------------------------------------------------------------------
	} else {
	#+++++++++++++++++++++ ONLY FOR HMM2 MODE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		$grpdetailsR = $grpid_statusR->{ $grpid}{"current_hmm_details"}; # use details from status file
		if (exists $grpid_sp_difscoresR->{$grpid}) { # some better matches
			#+ reorganise alt scores to be keyed by altid
			my %sp_difcores = %{ $grpid_sp_difscoresR->{$grpid} }; 
			foreach my $sp (keys %sp_difcores) {
				foreach my $altid (keys %{ $sp_difcores{$sp} } ) {
					$altid_sp{$altid} = $sp;
					$altid_scoredif{$altid} = $sp_difcores{$sp}{$altid};
				}
			}
			#-
			$altfaF = $faFsR->{"HMM1alt"};	
		}
	#-----------------------------------------------------------------------------------------------------------------------------------------
	}
	
	my $changed = 0;
	my $na = keys %altid_scoredif;
	if ($na > 0) {
		print "$na alt peps score better than existing members in $grpid. Starting optimisation\n";

		foreach my $altid (sort {$altid_scoredif{$b} <=> $altid_scoredif{$a}} keys %altid_scoredif) {
			my $sp = $altid_sp{$altid};
			print "\taltid $altid sp $sp scoredif ". $altid_scoredif{$altid} ."\n";
		}
		
		
		#--------- trial replacement of grp members
		#+ ++++++++ try any alt ids that scored better with most promising first
		foreach my $altid (sort {$altid_scoredif{$b} <=> $altid_scoredif{$a}} keys %altid_scoredif) {
			my $sp = $altid_sp{$altid};
			my ($isbetter, $resulttxt, $newgrpdetailsR) = HMMs::try_replacing_seq_in_group($grpid, $grpdetailsR, $sp, $altfaF, $altid);
			if ($isbetter) {
				$changed++;
				$grpdetailsR = $newgrpdetailsR;
			}
			print "\t$resulttxt\n";
		}
		# #- --------
	} else {
		print "\t(No alt sequences to try for $grpid)\n";
	}
	print "$steptodo complete $grpid. $changed improvements\n";

	# # Tried all alternatives. Optimisation finished. 
	$grpid_statusR->{ $grpid}{"current_hmm_details"} = $grpdetailsR;
	$grpid_statusR->{ $grpid}{"step_complete"} = $steptodo;
	Pipeline::update_fsteps_status($grpid_statusR, $fsstatF);
	$nb++ if ($changed>0);
}
print "\nEnd of step from ". $STEPTODO_STEPDONE{$steptodo} ." to $steptodo for $ng groups. $nb improved\n";

exit(0);
