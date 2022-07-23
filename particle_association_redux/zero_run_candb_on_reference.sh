#!/bin/bash
n=8
prefix=/workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces
dbdir=/workspace/data/zhongj/Transposase_Project/particle_lifestyle/db
pieces="/workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dd.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_de.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_df.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dg.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dh.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_di.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dj.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dk.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dl.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dm.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dn.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_do.fna /workspace/data/zhongj/Transposase_Project/particle_lifestyle/OM-RGC_v2_reference_pieces/p_dp.fna"
# pieces=$(ls $prefix/p_*s.fna) # opqrstuvwxyz
for p in $pieces; do
	sample=${p#"${prefix}/"}
	sample=${sample%%.fna}
	echo $sample
	inputfile=${prefix}/dbcan_output_old/${sample}/uniInput
	dir=${prefix}/dbcan_output/${sample}
	echo $dir "started at $(date)"
	run_dbcan $inputfile protein --out_dir $dir --dia_cpu $n --hmm_cpu $n --tf_cpu $n --stp_cpu $n --db_dir $dbdir
done
