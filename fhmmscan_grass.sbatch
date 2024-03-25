#!/bin/bash
## 
#SBATCH --job-name=fhmms_grass
#SBATCH --output=%j.fhmms_grass.out
#SBATCH --error=%j.fhmms_grass.err
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
module load BioPerl
module load MUSCLE
module load HMMER

pids=""
# get each grass sp from this file as query sp
input="grass_sp_status.txt"
while IFS= read -r line
do
	IFS=$'\t'
	tmp=($line)
	if [[ ${tmp[0]} != "#"* ]] && [[ ${tmp[0]} == *"_"* ]]; then
		echo "sbatch sp read: ${tmp[0]}"
		./hmmscan_allgm.pl ${tmp[0]} final_hmm &
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
		echo "All processes ok. Running script to analyse hits"
		./analyse-hmmscan-grasshits.pl
else
		echo "At least 1 process FAILed check error log"
fi

