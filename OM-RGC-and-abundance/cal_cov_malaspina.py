#!/usr/bin/env python3 
# Written by Jimmy Zhong (zhongj2@carleton.edu), Carleton '23 under Professor Rika Anderson
# date: Aug 18th 2022
import os
import glob
import pandas as pd
import csv
import statistics
# RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
# = total read base count in a MAG / ( geneLength/1000 * total read bases/1,000,000 )

base="/researchdrive/zhongj2/deep_ocean_bins"
base1=f"{base}/deep_metagenome_transposase_BLAST"

def get_target_gene_dict(blastp_file):
	gene_dict = {}
	with open(blastp_file) as file:
		blastp_res = list(csv.reader(file, delimiter="\t"))[1:]
	for l in blastp_res:
		if len(l) == 12:
			anvio_gene_caller_id, start_pos, end_pos = l[1], int(l[8]), int(l[9])
			if anvio_gene_caller_id not in gene_dict:
				gene_dict[anvio_gene_caller_id] = [abs(end_pos - start_pos)]
			else:
				# if there are mutliple hits to the same ORF, will take the mean gene length 
				gene_dict[anvio_gene_caller_id].append(abs(end_pos - start_pos))
		else:
			print("is this the last line?")
	for id in gene_dict.keys():
		len_list = gene_dict[id]
		median_gene_len = statistics.mean(len_list)
		gene_dict[id] = median_gene_len
	return gene_dict

def get_COG_gene_dict(cog_file, keyword, anvi_gene_call):
	defense_gene_id = set()
	with open(cog_file) as f:
		while True:
			line = f.readline()
			if not line: break
			l = line.split("\t")
			if keyword in l[3]:
				defense_gene_id.add(l[0])
	defense_dict = {}
	for l in anvi_gene_call:
		if l[0] in defense_gene_id:
			defense_dict[l[0]] = abs(int(l[3])-int(l[2]))
	return defense_dict

def get_CAZpep_dict(anvio_call_id, anvi_gene_call):
	gene_dict = {}
	gene_id = set(line.strip() for line in open(anvio_call_id))
	for l in anvi_gene_call:
		if l[0] in gene_id:
			gene_dict[l[0]] = abs(int(l[3])-int(l[2]))
	return gene_dict

def find_cov(gene_covs, gene_dict, total_bp):
	sum_cov = 0
	for line in gene_covs:
		l = line.split("\t")
		if len(l) < 2: continue
		anvio_gene_id, this_ORF_ave_cov = str(l[0]), float(l[1])
		if anvio_gene_id in gene_dict:
			sum_cov += this_ORF_ave_cov * gene_dict[anvio_gene_id]
	gene_prop = sum_cov / total_bp
	return sum_cov, gene_prop

def streamline_gene(gene_dict, new_col_name, df):
	md_list = df[["sample", "lower_filter_size", "Sequencing_depth_Gbp"]].values.tolist()
	gene_prop_list = []
	for l in md_list:
		with open(f"{base}/per-sample-mapping/{l[0]}_{l[1]}-GENE-COVERAGES.txt", "r") as f:
			gene_covs = f.readlines()[1:]
		total_bp = float(l[2]) * 1000000000
		gene_sum_cov, gene_prop = find_cov(gene_covs, gene_dict, total_bp)
		gene_prop_list.append(gene_prop)
	df[new_col_name] = gene_prop_list
	return df

anvi_gene_call = list(csv.reader(open(f"{base}/anvio-gene-call.txt"), delimiter='\t'))
anvi_category = f"{base}/all-anvi-cog-category.tsv"

trans_dict = get_target_gene_dict(f"{base1}/deep_anvi-genecall_transposase.blastp")
defense_dict = get_COG_gene_dict(anvi_category, "Defense mechanisms", anvi_gene_call)
signalT_dict = get_COG_gene_dict(anvi_category, "Signal transduction mechanisms", anvi_gene_call)
replication_dict = get_COG_gene_dict(anvi_category, "Replication, recombination and repair", anvi_gene_call)
toxin_dict = get_COG_gene_dict(f"{base}/all-anvi-cog-function.tsv", "toxin", anvi_gene_call)
all_CAZ_dict = get_CAZpep_dict(f"{base1}/all-CAZyme-norm.txt", anvi_gene_call)
sec_CAZ_dict = get_CAZpep_dict(f"{base1}/deep_secretory_CAZyme.txt", anvi_gene_call)
all_pep_dict = get_CAZpep_dict(f"{base1}/all-peptidase-norm.txt", anvi_gene_call)
sec_pep_dict = get_CAZpep_dict(f"{base1}/deep_secretory_peptidase.txt", anvi_gene_call)
amino_dict = get_COG_gene_dict(anvi_category, "Amino acid transport and metabolism", anvi_gene_call)
lipid_dict = get_COG_gene_dict(anvi_category, "Lipid transport and metabolism", anvi_gene_call)
coenzyme_dict = get_COG_gene_dict(anvi_category, "Coenzyme transport and metabolism", anvi_gene_call)
cwall_dict = get_COG_gene_dict(anvi_category, "Cell wall", anvi_gene_call)
energy_dict = get_COG_gene_dict(anvi_category, "Energy production and conversion", anvi_gene_call)

cmobility_dict = get_COG_gene_dict(anvi_category, "Cell motility", anvi_gene_call)
extstruct_dict = get_COG_gene_dict(anvi_category, "Extracellular structures", anvi_gene_call)
cytoskeleton_dict = get_COG_gene_dict(anvi_category, "Cytoskeleton", anvi_gene_call)
carbohydrate_dict = get_COG_gene_dict(anvi_category, "Carbohydrate transport and metabolism", anvi_gene_call)
inorgalIon_dict = get_COG_gene_dict(anvi_category, "Inorganic ion transport and metabolism", anvi_gene_call)
nucl_dict = get_COG_gene_dict(anvi_category, "Nucleotide transport and metabolism", anvi_gene_call)

all_dict = [trans_dict, defense_dict, signalT_dict, toxin_dict, replication_dict,
			all_CAZ_dict, sec_CAZ_dict, all_pep_dict, sec_pep_dict,
			amino_dict,lipid_dict,coenzyme_dict,cwall_dict,energy_dict]

str1 = "Transposase,Defense_mechanisms,Signal_transduction_mechanisms,ToxinAntiT"
str2 = "Replication_recombination_and_repair,CAZyme,secretory_CAZyme,peptidase,secretory_peptidase"
str3 = "Amino_acid_TM,Lipid_TM,Coenzyme_TM,Cell_wall,Energy_production_and_conversion"

cols = f"{str1},{str2},{str3}".split(",")
cols += ["Cell_motility", "Extracellular_struct", "Cytoskeleton", "Carbohydrate_TM", "Inorganic_ion_TM", "Nucleotide_TM"]
cols += ["Chromatin_structure_and_dynamics", "Cell_cycle_control", "Intracellular_trafficking", "Secondary_metabolites"]
cols_dna = ["DNA_"+ x for x in cols]

print("finish parsing dictionaries")

metadata = pd.read_csv(f"{base1}/sample_metadata.tsv", sep = "\t")
md_cols = "sample,station,lower_filter_size,Ocean_DNA,Ocean,Depth,Longitude,Latitude,Salinity_PSU,Temperature,Oxygen_DNA,WATERMASS,Basin,Sequencing_depth_Gbp".split(",")
metadata.columns = md_cols

for i in range(len(all_dict)):
	print(f"doing {cols_dna[i]}")
	metadata = streamline_gene(all_dict[i], cols_dna[i], metadata)

out_f = 'Malaspina-genes-coverage.csv'
metadata.to_csv(out_f, index=False)

# for things in the future:
# amino_dict = get_COG_gene_dict(anvi_category, "Amino acid transport and metabolism", anvi_gene_call)

chrotid_dict = get_COG_gene_dict(anvi_category, "Chromatin structure and dynamics", anvi_gene_call)
ccycle_dict = get_COG_gene_dict(anvi_category, "Cell cycle control", anvi_gene_call)
intracellular_dict = get_COG_gene_dict(anvi_category, "Intracellular trafficking", anvi_gene_call)
sec_metabolite_dict = get_COG_gene_dict(anvi_category, "Secondary metabolites", anvi_gene_call)

out_f = 'Malaspina-genes-coverage.csv'
metadata = pd.read_csv(out_f)

# dicts =[chrotid_dict, ccycle_dict, intracellular_dict, sec_metabolite_dict]
# w_cols = ["Chromatin_structure_and_dynamics", "Cell_cycle_control", "Intracellular_trafficking", "Secondary_metabolites"]
# words = ["DNA_"+ x for x in w_cols]

# for i in range(len(dicts)):
# 	print("doing " + words[i])
# 	metadata = streamline_gene(dicts[i], words[i], metadata)

# metadata = streamline_gene(extstruct_dict, "DNA_CAZyme", metadata)
# metadata = streamline_gene(cytoskeleton_dict, "DNA_peptidase", metadata)

metadata.to_csv(out_f, index=False)
