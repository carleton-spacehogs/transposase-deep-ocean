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

samples=$(ls /researchdrive/zhongj2/tara_contig_fasta/Tara*/*.bam) # only 0.22-1.6

for s in $samples; do
bamfile=$(echo $s | awk -F'/' '{print $(NF)}')
outf=${s%%${bamfile}}contig_cov.txt
samtools idxstats $s > $outf
newf=$(echo $bamfile | awk -F'_' '{print $2}')_contig_cov.txt
cp $outf $newf
echo $newf
done

