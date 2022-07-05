#!/bin/bash

# working directory: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/run_dbcan/tools/signalp-5.0b/signalp

aminof=$(ls ../../../../OM-RGC_v2_reference_pieces/dbcan_output_old/p_*/only_CAZenzyme.faa)
for f in $aminof; do
sample=$(echo $f | awk -F'/' '{print $7}')
./signalp -org 'gram-' -fasta $f -prefix ${sample}_gramNegative
./signalp -org 'gram+' -fasta $f -prefix ${sample}_gramPositive
done
