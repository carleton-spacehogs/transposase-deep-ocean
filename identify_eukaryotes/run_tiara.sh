#!/bin/bash

source activate tiara_euk

base="/researchdrive/zhongj2/deep_ocean_bins/per-assembly-sample-mapping"
size=02 # 0.8 - 5 mu-meter; 0.2 - 0.8 mu-meter
samples=$(ls ${base}/per-assembly-per-sample-mapping${size}/*/*_min1000.fasta)

# /researchdrive/zhongj2/deep_ocean_bins/per-assembly-sample-mapping/per-assembly-per-sample-mapping08/SRR3960580/SRR3960580_min1000.fasta
for s in $samples; do
outf=$(echo $s | awk -F'/' '{print $(NF)}')
outf=${outf%%.fasta}_size${size}_euk.txt
tiara -i $s -o $outf -t 10
done

size=3
samples=$(ls /researchdrive/zhongj2/tara_contig_fasta/Tara*/*${size}.fasta)
for s in $samples; do
news=${s%%.fasta}_min1000.fasta
echo $news
# ~/seqtk/seqtk seq -L 1000 $s > $news
outf=$(echo $news | awk -F'/' '{print $(NF)}')
outf=${outf%%.fasta}_euk.txt
tiara -i $s -o $outf -t 20
done

