#!/bin/bash

# working directory: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/run_dbcan/tools/signalp-5.0b/bin
# to run: ../../../../transposase-deep-ocean/particle-association-redux/second_signalp_find_secretory.sh

# aminof=$(ls ../../../../OM-RGC_v2_reference_pieces/dbcan_output_old/p_*/only_CAZenzyme.faa)
# for f in $aminof; do
# sample=$(echo $f | awk -F'/' '{print $7}')
# ./signalp -org 'gram-' -fasta $f -prefix ${sample}_CAZenzyme_gramNegative
# ./signalp -org 'gram+' -fasta $f -prefix ${sample}_CAZenzyme_gramPositive
# done

aminofp=$(ls ../../../../OM-RGC_v2_reference_pieces/dbcan_output_old/p_d*/only_peptidase.faa)
for f in $aminofp; do
sample=$(echo $f | awk -F'/' '{print $7}')
echo "I am on ${f}"
./signalp -org 'gram-' -fasta $f -prefix ${sample}_gramNegative
./signalp -org 'gram+' -fasta $f -prefix ${sample}_gramPositive
done

