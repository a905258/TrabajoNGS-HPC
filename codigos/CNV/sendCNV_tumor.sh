#!/bin/bash

files=("SRR645209" "SRR645218" "SRR645222" "SRR645224" "SRR645235" "SRR645241" \
"SRR645243" "SRR645247" "SRR645249" "SRR645259" "SRR645272" "SRR645274" \
"SRR645283" "SRR645289")
nums=${#files[@]}
num=$((nums-1))
echo $num
sbatch --array=0-$num%10 /scratch/a905383/CNV/CNV_tumor.sbs

