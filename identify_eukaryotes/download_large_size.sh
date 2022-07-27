#!/bin/bash

depth_samples=$(grep "5 - 20" Carradec-et-al-AGlobalOceanAtlasOfEukaryoticGenes-SuppTable6.csv | grep ${1} | shuf -n 20 | sed 's/ //g')
> ${1}_5-20_sampled.txt
pre="ascp -QT -l 300m -P 33001 -i /Accounts/zhongj2/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq"
root="/researchdrive/zhongj2/MAG_pnps_redux"
depth_folder=${root}/${1}_5-20/raw_reads

mkdir $depth_folder

for line in $depth_samples; do
echo $line >> ${1}_5-20_sampled.txt
ERRs=$(echo $line | awk -F"," '{print $12}' | sed s/\"//g)
echo $ERRs
for ERR in $ERRs; do 
ERR_pre=$(echo $ERR | cut -c 1-6)
${pre}/${ERR_pre}/00${ERR: -1}/${ERR}/${ERR}_1.fastq.gz $depth_folder &
${pre}/${ERR_pre}/00${ERR: -1}/${ERR}/${ERR}_2.fastq.gz $depth_folder
done
done
