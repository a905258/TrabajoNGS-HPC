#!/bin/bash
###################################################################################
### CNV
###################################################################################
#
#

#SBATCH --job-name=cnv%a
#SBATCH --account=biomedicine
#SBATCH --partition=biomedicine
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --time=36:30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=a905383@alumni.unav.es
#SBATCH -o /scratch/a905383/CNV/logs/names%a.log


### Load programs here

file=($(cd /scratch/a905383/CNV/counts && basename --suffix=.denoisedCR.tsv -- *.denoisedCR.tsv))
filename=${file[$SLURM_ARRAY_TASK_ID]}

bash /scratch/a905383/CNV/igv_files/formater.sbs "/scratch/a905383/CNV/counts/"$filename".denoisedCR.tsv"

