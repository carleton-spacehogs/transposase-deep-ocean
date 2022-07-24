
root=/workspace/data/zhongj/Transposase_Project/particle_lifestyle
/researchdrive/zhongj2/psortb -a --seq ${root}/bins/ARS/TOBG_ARS-1018.fna --outdir .

/Accounts/zhongj2/psortb/bin/psort -a $root/bins/ARS/TOBG_ARS-1018/only_CAZenzyme.faa --output .

conda install -c bioconda perl-bioperl
perl Makefile.PL

/Accounts/zhongj2/bio-tools-psort-all/blast-2.2.17/bin/blastall

pfscan: ~/miniconda3/bin/pfam_scan.pl

/Accounts/zhongj2/hmmer/bin
~/miniconda3/bin/squid

PSORT_ROOT
PSORT_PFTOOLS
BLASTDIR
