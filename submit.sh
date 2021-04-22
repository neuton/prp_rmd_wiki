#!/bin/bash
function sbatch_id { sbatch $* | cut -d" " -f4; }

unfold=`sbatch_id unfold.sh`
echo $unfold

for ((i=1; i<=10; i++)); do
	sbatch_id --dependency=afterok:$unfold equilibrate.sh $i
done