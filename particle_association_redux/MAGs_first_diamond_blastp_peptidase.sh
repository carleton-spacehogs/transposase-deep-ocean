#!/bin/bash

source activate # base environment has diamond

# working directory: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux
# to run: ./MAGs_first_diamond_blastp_peptidase.sh

peptidase_db="../../pepunit.lib"
# ocean="ARS"
# ocean="CPC"
# ocean="IN"
# ocean="NP"
# ocean="RS"
# ocean="EAC"
# ocean="MED"
# ocean="SP"
# all_MAGs=$(ls -d ../../bins/deep/mp-deep_mag-*/)
# ocean="NAT"
ocean="SAT"

all_MAGs=$(ls -d ../../bins/${ocean}/*/)

for MAG in $all_MAGs; do
diamond_db=${MAG}MAG_aa_seq
blastp_out=${MAG}peptidase_diamond_unique.blastp
diamond makedb --in ${MAG}uniInput -d $diamond_db -p 5
# https://unix.stackexchange.com/questions/29961/get-lines-with-maximum-values-in-the-column-using-awk-uniq-and-sort
diamond blastp -d $diamond_db -q $peptidase_db --threads 5 --evalue 0.0000000001 | sort -k11 | sort -u -k2,2 > $blastp_out
wc -l $blastp_out
done

