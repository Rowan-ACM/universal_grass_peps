#!/bin/bash
# 
#SBATCH --job-name=optHMM1
#SBATCH --output=%j.optHMM1.out
#SBATCH --error=%j.optHMM1.err
#SBATCH --time=7-0
#SBATCH --partition=standard
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16384
#SBATCH --cpus-per-task=8
#SBATCH --spread-job
#SBATCH --get-user-env=L
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rowan.mitchell@rothamsted.ac.uk
#+
# Below here is the script itself
#-
module load BioPerl
module load MUSCLE
module load HMMER

./split-pipeline-statfiles.pl

# Run 12 concurrent processes
pids=""
for i in {1..12}; do
        ./optimise_hmms.pl HMM1 $i &
        # store PID of process
        pids+=" $!"
done

# Wait for all processes to finish
i=0
for p in $pids; do
		((i++))
        if wait $p; then
                echo "Process # $i ID $p success"
        else
                echo "Process # $i ID $p fail"
        fi
done

./join-pipeline-statfiles.pl HMM1
