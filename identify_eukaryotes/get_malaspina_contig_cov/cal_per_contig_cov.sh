#!/bin/bash

base="/researchdrive/zhongj2/deep_ocean_bins/per-assembly-sample-mapping"
samples=$(ls ${base}/per-assembly-per-sample-mapping*/*/*_sorted.bam)

for s in $samples; do
bamfile=$(echo $s | awk -F'/' '{print $(NF)}')
outf=${s%%${bamfile}}contig_cov.txt
samtools idxstats $s > $outf
cp $outf ${bamfile%%_sorted.bam}_contig_cov.txt
echo $s $outf
done

