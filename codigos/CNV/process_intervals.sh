source /scratch/arubio/bin/moduleloadNGS

#gatk PreprocessIntervals -R /scratch/arubio/DNAseq/Auxfiles/dummy.fa \
#-L 1 \
#--bin-length 1000 --interval-merging-rule OVERLAPPING_ONLY \
#--padding 0 -O /scratch/a905383/CNV/preprocessed_intervals1.interval_list

gatk FilterIntervals -L /scratch/a905383/CNV/preprocessed_intervals10000.interval_list \
 --interval-merging-rule OVERLAPPING_ONLY \
 -XL blacklist.interval_list \
 -I /scratch/a905383/CNV/counts/SRR645241.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645209.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645247.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645222.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645272.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645218.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645243.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645283.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645274.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645224.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645249.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645289.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645235.counts.hdf5 \
 -I /scratch/a905383/CNV/counts/SRR645259.counts.hdf5 \
 --annotated-intervals /scratch/a905383/CNV/annotated_intervals.tsv \
 --minimum-gc-content 0.54 \
 -O /scratch/a905383/CNV/filtered.interval_list





