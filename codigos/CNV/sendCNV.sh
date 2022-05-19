#!/bin/bash

files=("SRR645210" "SRR645219" "SRR645223" "SRR645225" "SRR645236" "SRR645242" \
"SRR645244" "SRR645248" "SRR645250" "SRR645260" "SRR645273" "SRR645275" \
"SRR645284" "SRR645288")

nums=${#files[@]}
num=$((nums-1))
echo $num
sbatch --array=0-$num%10 /scratch/a905383/CNV/CNV.sbs

