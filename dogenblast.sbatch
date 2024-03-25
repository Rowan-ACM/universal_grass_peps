#!/bin/bash
# 
#SBATCH --job-name=dogenblast
#SBATCH --output=%j.dogenblast.out
#SBATCH --error=%j.dogenblast.err
#SBATCH --time=3-12
#SBATCH --partition=standard
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=16384
#SBATCH --cpus-per-task=8
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
        ./do-genblast.pl $i &
        # store PID of process
        pids+=" $!"
done

# Wait for all processes to finish
allok=1
i=0
for p in $pids; do
		((i++))
        if wait $p; then
                echo "Process # $i ID $p success"
        else
				allok=0
               echo "Process # $i ID $p fail"
        fi
done

# re-join stat files regardless of whether all processes succeeded
./join-pipeline-statfiles.pl finalHMM
