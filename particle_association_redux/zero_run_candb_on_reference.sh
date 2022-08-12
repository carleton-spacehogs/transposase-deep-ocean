#!/bin/bash
n=8
prefix=/workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces
dbdir=/workspace/data/zhongj/Transposase_Project/particle_lifestyle/db
pieces=$(echo {o..o})
# pieces=$(ls $prefix/p_*s.fna) # opqrstuvwxyz
for p in $pieces; do
	sample=p_c${p}
	# sample=${p#"${prefix}/"}
	# sample=${sample%%.fna}
	echo $sample
	inputfile=${prefix}/dbcan_output_old/${sample}/uniInput
	dir=${prefix}/dbcan_output/${sample}
	echo $dir "started at $(date)"
	run_dbcan $inputfile protein --out_dir $dir --dia_cpu $n --hmm_cpu $n --tf_cpu $n --stp_cpu $n --db_dir $dbdir
done
