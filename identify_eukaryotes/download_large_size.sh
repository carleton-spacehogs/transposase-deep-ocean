#!/bin/bash

low=20
high=180
depth=${1}

depth_samples=$(grep "$low - $high" Carradec-et-al-AGlobalOceanAtlasOfEukaryoticGenes-SuppTable6.csv | grep $depth | shuf -n 2 | sed 's/ //g')
pre="ascp -QT -l 300m -P 33001 -i /Accounts/zhongj2/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq"
depth_folder=/researchdrive/SpaceHogs_shared/jimmy_stuff/DNA_raw_reads/${depth}_${low}_${high}_reads

mkdir $depth_folder

for line in $depth_samples; do
echo $line >> ${depth}_${low}_${high}_sampled.txt
ERRs=$(echo $line | awk -F"," '{print $12}' | sed s/\"//g)
echo $ERRs
for ERR in $ERRs; do 
ERR_pre=$(echo $ERR | cut -c 1-6)
${pre}/${ERR_pre}/00${ERR: -1}/${ERR}/${ERR}_1.fastq.gz $depth_folder &
${pre}/${ERR_pre}/00${ERR: -1}/${ERR}/${ERR}_2.fastq.gz $depth_folder
done
done
