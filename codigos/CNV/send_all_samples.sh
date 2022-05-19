#!/bin/bash

file=($(ls /scratch/a905383/CNV/counts/*.denoisedCR.tsv | wc -l))
num=$((file-1))

echo $num
sbatch --array=0-$num%10 /scratch/a905383/CNV/igv_files/all_samples.sh


