#!/bin/bash

file=($(ls /scratch/a905383/CNV/raw_data/*_aligned.bam | wc -l))
num=$((file-1))

echo $num
sbatch --array=0-$num%10 /scratch/a905383/CNV/secondStep.sbs

