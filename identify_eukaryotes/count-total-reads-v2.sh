#!/bin/bash

# /researchdrive/zhongj2/MAG_pnps_redux/IN/raw_reads
# /researchdrive/zhongj2/MAG_pnps_redux/SRF_5-20/raw_reads/
# /workspace/data/zhongj/Transposase_Project/other_DNA_to_keep
ocean=$1

prefix="/researchdrive/zhongj2/MAG_pnps_redux/$ocean/raw_reads/"
suffix="_1.fastq.gz"
samples=$(ls ${prefix}*${suffix})
for sample in $samples; do
if [ -f "$sample" ]; then
	s_sample=${sample#"$prefix"}
	s_sample=${s_sample%"$suffix"}
	echo counting sample: $s_sample
	num_lines=$(zcat $sample|wc -l)
	echo $num_lines
	num_reads=$(echo "$num_lines/2" | bc) # /4 but times 2 (forward and backward reads)
	echo $s_sample,$num_reads >> count-${ocean}-reads.out 
else
	echo cannot find sample: $sample
fi
done