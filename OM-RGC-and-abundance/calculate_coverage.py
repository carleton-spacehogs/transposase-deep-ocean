import pandas as pd
import numpy as np
import csv
reference = list(csv.reader(open("OM-RGC_v2.tsv"), delimiter='\t'))

def get_id(blast_file):
	# "OM-RGC_v2_Reference_Pieces/Biofilm_BLAST_Result/merged_Biofilm_hits.txt"
	biofilm_blast = list(csv.reader(open(blast_file), delimiter='\t'))
	biofilm_id = set()
	for line in biofilm_blast:
		OMRGC_ID = line[1].split("_")[2] # 'TARA_Y200000002_OM-RGC.v2.004389382_COG3191_Bacteria'
		biofilm_id.add(OMRGC_ID)
	return biofilm_id

def get_gene_specific_ref(gene_id_set, gene_name):
	out = []
	gene_length_dict = {}
	for line in reference:
		if line[1] in gene_id_set:
			gene_id, COG, gene_length = line[1], line[3], len(line[13])
			out.append(f"{gene_id},{COG},{gene_length}")
			gene_length_dict[gene_id] = gene_length
	with open(f"{gene_name}_OMRGC_ID_length.csv", 'w') as f:f.write("\n".join(out))
	return gene_length_dict

def get_gene_specific_ref_through_COG(search_string):
	cog_out = set()
	cogfile = open('/usr/local/data/Python_scripts/COG_numbers_categories_annotations.txt', 'r')
	for line in cogfile:
		columns = line.split('\t')
		cog_number, cog_category = columns[0], columns[2]
		if cog_category == search_string: 
			cog_out.add(cog_number)
	cogfile.close()
	out = []
	gene_length_dict = {}
	for line in reference:
		if line[3] in cog_out:
			gene_id, COG, gene_length = line[1], line[3], len(line[13])
			out.append(f"{gene_id}, {COG}, {gene_length}")
			gene_length_dict[gene_id] = gene_length
	with open(f"{search_string}_OMRGC_ID_length.csv", 'w') as f:f.write("\n".join(out))
	return gene_length_dict

def subset_coverage_rows(gene_length_dict, gene_name, coverage_file, DNA_or_RNA):
	counter = 0
	out = []
	# coverage_file = 'OM-RGC_v2_gene_profile_metaG.tsv', or '../Tara_metatranscriptomic/OM-RGC_v2_gene_profile_metaT.tsv'
	with open(coverage_file) as coverage:
		out.append(coverage.readline().rstrip() + "\tgene_length") # firstline 
		for line in coverage:
			counter += 1
			gene_id = line[0:19]
			if gene_id in gene_length_dict:
				gene_length=f"\t{str(gene_length_dict[gene_id])}"
				newline = line.strip() + gene_length
				out.append(newline)
				if (counter % 100000 == 0): 
					print("I am on line: " + str(counter))
	with open(f"{gene_name}_{DNA_or_RNA}_coverage_subset.csv", 'w') as f:f.write("\n".join(out))

# abandoned because too slow and took up too much memory
# def cal_sum_num_reads_mapped(gene_name, DNA_or_RNA):
# 	subset = pd.read_csv(f"{gene_name}_{DNA_or_RNA}_coverage_subset.csv", sep = "\t")
# 	print(subset)
# 	length_normalized = subset.multiply(subset["gene_length"], axis="index")
# 	length_normalized = length_normalized.drop('gene_length', axis=1)
# 	sum_length = length_normalized.sum(axis = 0, skipna = True)
# 	print(f"finish calculating the sum_length for {gene_name}")
# 	return sum_length

import operator
def cal_sum_num_reads_mapped(gene_name, DNA_or_RNA):
	counter = 0
	subset = list(csv.reader(open(f"{gene_name}_{DNA_or_RNA}_coverage_subset.csv") , delimiter='\t'))
	col_num = 187 if DNA_or_RNA == "RNA" else 180
	sum_row = [0]*col_num
	for row in subset[1:]:
		gene_len, abundance = int(row[-1]), row[1:-1]
		float_abundance = [float(item) for item in abundance]
		counter += 1 
		if (counter % 10000 == 0): 
			print("I am on line: " + str(counter))
		result = list(map(lambda x : x * gene_len, float_abundance)) # sum_row = list(map(lambda x : x + 5, sum_row))
		sum_row = list(map(operator.add, result, sum_row)) # list(map(operator.add, sum_row, sum_row))
	return_df = [subset[0][1:-1], sum_row]
	return pd.DataFrame(return_df).transpose()


def cal_prop_reads_mapped(sum_length, gene_name, DNA_or_RNA):
	sum_len = sum_length.set_axis(["sample", "sum_bases"], axis=1)
	# sum_len = sum_length.reset_index().set_axis(["sample", "sum_bases"], axis=1) # for the old version cal_sum_num_reads_mapped()
	normalization_factor = pd.read_csv(f"sample_readCounts_readLength.tsv", sep = "\t")
	if DNA_or_RNA == "RNA":
		normalization_factor = pd.read_csv(f"sample_readCounts_readLength_RNA.tsv", sep = "\t")
	sum_len = sum_len.merge(normalization_factor, on="sample", how="inner")
	print(sum_len)
	sum_len['read_counts'] = sum_len['sum_bases']/sum_len['Avg_read_len']
	sum_len['prop_reads_mapped_back_genes'] = sum_len['read_counts']/sum_len['HQ_reads']
	sum_len.to_csv(f'prop_{DNA_or_RNA}_reads_mapped_back_{gene_name}.csv', sep=',', index=False)
	return sum_len

# biofilm_id_set = get_id("OM-RGC_v2_Reference_Pieces/Biofilm_BLAST_Result/merged_Biofilm_hits.txt")
biofilm_id_set = get_id("OM-RGC_v2_Reference_Pieces/Biofilm_BLAST_Result/merge_Biofilm_hits_no_BapA.txt")
biofilm_length_dict = get_gene_specific_ref(biofilm_id_set, "biofilm")
trans_id_set = get_id("OM-RGC_v2_Reference_Pieces/Transposase_BLAST_Result/merged_Transposase_hits.txt")
trans_length_dict = get_gene_specific_ref(trans_id_set, "transposase")
defense_length_dict = get_gene_specific_ref_through_COG("Defense_mechanisms")

# read existing dicts
def read_dict(gene_name):
	id_length_dict = {}
	dict_file = list(csv.reader(open(f"{gene_name}_OMRGC_ID_length.csv"), delimiter=','))
	for line in dict_file:
		id_length_dict[line[0].strip()] = line[2].strip() # f"{gene_id}, {COG}, {gene_length}"
	return id_length_dict

biofilm_length_dict, trans_length_dict, defense_length_dict = read_dict("biofilm"), read_dict("transposase"), read_dict("Defense_mechanisms")

print("got all the dicts")
# DNA
subset_coverage_rows(biofilm_length_dict, "biofilm", "OM-RGC_v2_gene_profile_metaG.tsv", "DNA")
subset_coverage_rows(trans_length_dict, "transposase", "OM-RGC_v2_gene_profile_metaG.tsv", "DNA")
subset_coverage_rows(defense_length_dict, "Defense_mechanisms", "OM-RGC_v2_gene_profile_metaG.tsv", "DNA")
biofilm_sum = cal_sum_num_reads_mapped("biofilm", "DNA")
biofilm_prop = cal_prop_reads_mapped(biofilm_sum, "biofilm", "DNA")
transposase_sum = cal_sum_num_reads_mapped("transposase", "DNA")
transposase_prop = cal_prop_reads_mapped(transposase_sum, "transposase", "DNA")
defense_sum = cal_sum_num_reads_mapped("Defense_mechanisms", "DNA")
defense_prop = cal_prop_reads_mapped(defense_sum, "Defense_mechanisms", "DNA")

print("got all the DNA abundance, calculating the ones for RNA")
# RNA
subset_coverage_rows(biofilm_length_dict, "biofilm", "../Tara_metatranscriptomic/OM-RGC_v2_gene_profile_metaT.tsv", "RNA")
subset_coverage_rows(trans_length_dict, "transposase", "../Tara_metatranscriptomic/OM-RGC_v2_gene_profile_metaT.tsv", "RNA")
subset_coverage_rows(defense_length_dict, "Defense_mechanisms", "../Tara_metatranscriptomic/OM-RGC_v2_gene_profile_metaT.tsv", "RNA")
biofilm_sum = cal_sum_num_reads_mapped("biofilm", "RNA")
biofilm_prop = cal_prop_reads_mapped(biofilm_sum, "biofilm", "RNA")
transposase_sum = cal_sum_num_reads_mapped("transposase", "RNA")
transposase_prop = cal_prop_reads_mapped(transposase_sum, "transposase", "RNA")
defense_sum = cal_sum_num_reads_mapped("Defense_mechanisms", "RNA")
defense_prop = cal_prop_reads_mapped(defense_sum, "Defense_mechanisms", "RNA")

