#!/bin/bash

# /researchdrive/zhongj2/MAG_pnps_redux/IN/raw_reads
# /researchdrive/zhongj2/MAG_pnps_redux/SRF_5-20/raw_reads/
# /workspace/data/zhongj/Transposase_Project/other_DNA_to_keep
search_dir=$1

while read sample; do
	sample_path=${search_dir}/${sample}_1.fastq.gz
	if [ -f "$sample_path" ]; then
		echo counting sample: $sample
		num_lines=$(zcat $sample_path|wc -l)
		echo $num_lines
		num_reads=$(echo "$num_lines/2" | bc) # /4 but times 2 (forward and backward reads)
		echo $sample,$num_reads >> count-reads3.out 
	else
		echo cannot find sample: $sample_path
	fi
done < all_samples.txt
