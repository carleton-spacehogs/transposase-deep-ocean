#!/bin/bash
root=/workspace/data/zhongj/Transposase_Project/particle_lifestyle
demingr=/researchdrive/zhongj2/archaea-signalp

archMAGs=$(cat archaea-MAGs.txt)

for MAG in $archMAGs; do
	if [[ "$MAG" == *"deep_MAG_"* ]]; then
		MAG_folder=$(echo $MAG | sed 's/deep_MAG_/mp-deep_mag-/g')
		ocean="deep"
	else
		MAG_folder=TOBG_$MAG
		ocean=$(echo $MAG | awk -F'-' '{print $1}')
	fi
	echo $MAG_folder
	MAG_folder2=$root/bins/$ocean/$MAG_folder
	res_dir=$demingr/res-peptidase/$MAG
	mkdir $res_dir
	# $demingr/psortb -a --seq $MAG_folder2/only_CAZenzyme.faa --verbose --output terse --outdir $res_dir
	$demingr/psortb -a --seq $MAG_folder2/only_diamond_peptidase.faa --verbose --output terse --outdir $res_dir
	# cp $res_dir/*_psortb_archaea.txt $MAG_folder2/signal-CAZyme_psortb_archaea.txt
	cp $res_dir/*_psortb_archaea.txt $MAG_folder2/signal-peptidase_psortb_archaea.txt
done

# /Accounts/zhongj2/psortb/bin/psort -a $root/bins/ARS/TOBG_ARS-1018/only_CAZenzyme.faa --output .
# conda install -c bioconda perl-bioperl
# perl Makefile.PL
# /Accounts/zhongj2/bio-tools-psort-all/blast-2.2.17/bin/blastall
# pfscan: ~/miniconda3/bin/pfam_scan.pl
# /Accounts/zhongj2/hmmer/bin
# ~/miniconda3/bin/squid

