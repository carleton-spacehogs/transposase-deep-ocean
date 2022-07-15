#!/bin/bash

source activate # base environment has diamond

# working directory: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/toxin-db

toxin_db="all_toxin.fas"
oceans="IN NP SP NAT SAT deep MED EAC RS CPC ARS"

for ocean in $oceans; do
	all_MAGs=$(ls -d ../../bins/${ocean}/*/)
	for MAG in $all_MAGs; do
		diamond_db=${MAG}MAG_aa_seq
		blastp_out=${MAG}toxin_diamond_unique.blastp
		# already made: diamond makedb --in ${MAG}uniInput -d $diamond_db -p 5
		diamond blastp -d $diamond_db -q $toxin_db --more-sensitive --threads 5 --evalue 0.0000000001 | sort -k11 | sort -u -k2,2 > $blastp_out
		wc -l $blastp_out
	done
done

