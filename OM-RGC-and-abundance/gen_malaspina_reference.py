#!/usr/bin/env python3 
# Written by Jimmy Zhong (zhongj2@carleton.edu), Carleton '23 under Professor Rika Anderson
# date: Aug 28th 2022
import glob
import pandas as pd
import pysam
# output: the number of reads mapped to each ORFs
# note: malaspina read length: 150

base="/researchdrive/zhongj2/deep_ocean_bins"
base1=f"{base}/OM-RGC-ref"
oceans = ["SP", "SAT", "IN", "NAT", "NP"]
ocean = oceans[4] # 0, 1, 2, 3

## prep step: parse the reference into something useful
'''
ref = f"{base1}/M-geneDB.v1_1115269-seqs_full-header-annotation.txt"
ref_out = []
i = 0
for l in open(ref).readlines():
	if i % 10000 == 0:
		print(i)
	i += 1
	id = l[0:21] # M-geneDB.v1_001115261
	taxon, COG_num = "unTax", "unCOG"
	if "COG" in l:
		COG_num = l.split("_COG")[1]
		COG_num = COG_num.split("_")[0]
		COG_num = "COG" + COG_num # mutliple COG situation exists
		# i.e. _COG1148,COG1149_PF02662,PF12838_
	if "_k_" in l:
		taxon = l.split("_k_")[-1].strip()
	ref_out.append([id, COG_num, taxon])

ref_df = pd.DataFrame(ref_out)
ref_df.columns = ["id", "COG", "taxonomy"]
ref_df.to_csv("geneDB-annotation.csv", index=False)
'''
def as_db(num): # M-geneDB.v1_001115269
	return f"M-geneDB.v1_{str(num).zfill(9)}"

def count_iterable(i):
	return sum(1 for e in i)

def get_abun_list(bam_file):
	samfile = pysam.AlignmentFile(bam_file, "rb")
	sample_cov = []
	for id in geneIDs:
		if id % 10000 == 0:
			print(id)
		DB_id = as_db(id)
		num_reads = count_iterable(samfile.fetch(DB_id))
		sample_cov.append([DB_id, num_reads])
	sample_df = pd.DataFrame(sample_cov)
	sample_id = bam_file.split("/")[-1].replace("_sorted.bam", "")
	sample_df.columns = ["id", sample_id]
	return sample_df

ref_df = pd.read_csv(f"{base1}/geneDB-annotation.csv")
geneIDs = range(1, 1115269 + 1)
all_bams = glob.glob(f"{base1}/{ocean}/*_sorted.bam")

for bam_f in all_bams:
	sample_df = get_abun_list(bam_f)
	ref_df = ref_df.merge(sample_df, on = "id")
	print(ref_df)

ref_df.to_csv(f'{base1}/{ocean}/ref-abundance.csv', index=False)

# ref_df.iloc[:, 3:10].sum()*150/1000000000

''' merging the abundance matrix
start_df = pd.read_csv(f'{base1}/{oceans[0]}/ref-abundance.csv')
for ocean in oceans[1:]:
	cur_df = pd.read_csv(f'{base1}/{ocean}/ref-abundance.csv')
	start_df = start_df.merge(cur_df, on = ["id", "COG", "taxonomy"])
	print("merged: " + ocean)
start_df.to_csv(f'deep-ref-abundance.csv', index=False)
'''
df = pd.read_csv("deep-ref-abundance.csv")
reads_cov = df.iloc[:, 3:].sum()*150/1000000000
metadata = pd.read_csv(f"{base}/deep_metagenome_transposase_BLAST/sample_metadata.tsv", sep = "\t")

md = metadata[["ENA_Run_Accession_ID", "Sequencing_depth(Gbp)"]]
reads_cov = reads_cov.reset_index()
reads_cov.columns = ["ENA_Run_Accession_ID", "Mapped"]

res = md.merge(reads_cov, on = "ENA_Run_Accession_ID")
res['percent_mapped'] = res.Mapped / res["Sequencing_depth(Gbp)"]
