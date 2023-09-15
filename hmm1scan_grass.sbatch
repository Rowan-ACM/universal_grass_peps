#!/bin/bash
## 
#SBATCH --job-name=hmm1s_grass
#SBATCH --output=%j.hmm1s_grass.out
#SBATCH --error=%j.hmm1s_grass.err
#SBATCH --time=7-0
#SBATCH --partition=standard
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=32768
#SBATCH --cpus-per-task=8
#SBATCH --spread-job
#SBATCH --get-user-env=L
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rowan.mitchell@rothamsted.ac.uk
#+
# Below here is the script itself
#-
module purge
module load BioPerl
module load MUSCLE
module load HMMER

# make db
./hmmscan_allgm.pl none HMM1

pids=""
# get each grass sp from this file as query sp
input="grass_sp_status.txt"
while IFS= read -r line
do
	IFS=$'\t'
	tmp=($line)
	if [[ ${tmp[0]} != "#"* ]] && [[ ${tmp[0]} == *"_"* ]]; then
		echo "sbatch sp read: ${tmp[0]}"
		./hmmscan_allgm.pl ${tmp[0]} HMM1 &
        # store PID of process
        pids+=" $!"
	fi
done < "$input"

# Wait for all processes to finish and check all success
allok=1
IFS=$' '
for p in $pids; do
        if wait $p; then
                echo "Process $p success"
        else
				allok=0
                echo "Process $p fail"
        fi
done
if ((allok)); then
		echo "All processes ok. Running script to write fasta of hits"
		./write_hmm1scan_hits.pl
else
		echo "At least 1 process FAILed check error log"
fi

