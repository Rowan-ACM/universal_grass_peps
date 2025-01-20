The universal_grass_peps database doi.org/10.23637/rothamsted.98ywz is mainly for biologists to use and the scripts located here are not currently designed for easy re-use. Nevertheless with some changes it should be possible to run them on a Linux cluster.

You will need following modules available to install: BioPerl BLAST+ MUSCLE HMMER
And some way of running the genBlastG application.
 
Download fasta format files for all peps from gene models for all grass genomes, complete genome sequence files and genome annotation files in gff3 format. Also download ortholog predictions for rice and maize to each of the other grasses as separate tab-delimited files from Ensembl Plant Biomart. 

For later specificty steps download fasta format files for all peps from gene models for commelinid, monocot and non-monocot plant genomes to separate directory.

Change the hard coded directories in Pipeline.pm to your directory structure. 
genBlastG is tricky to run as it requires all files to be local and has some bugs that the scripts should correct.

Some quick scripts can be run from command line, most are run under slurm from sbatch files. The order of use and what they do is described in scripts.tsv

Rowan Mitchell, Jan 2025
