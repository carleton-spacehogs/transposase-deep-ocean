#!/bin/bash

source activate tiara_euk # base environment has diamond

# working directory: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/identify_eukaryotes
# to run: ./MAGs_find_euk_contigs.sh

# for ocean in SP SAT NAT deep; do
for ocean in ARS CPC IN NP RS EAC MED; do
	echo "currently working on $ocean ..." > logs/${ocean}.log
	bins=$(ls ../../bins/${ocean}/TOBG_${ocean}*.fna)
	for MAG in $bins; do
		MAG_name=$(echo $MAG | awk -F"/" '{print $5}')
		outf="../../bins/${ocean}/${MAG_name%%.fna}/tiara_euk.txt"
		echo tiara -i $MAG -o $outf -t 10 >> logs/${ocean}.log
		~/.local/bin/tiara -i $MAG -o $outf -t 10
	done
done

