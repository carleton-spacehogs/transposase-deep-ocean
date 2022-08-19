import pandas as pd
import numpy as np
import csv
import operator

base1="/researchdrive/zhongj2/BLAST_tara_OM-RGC_v2"
base2="/researchdrive/zhongj2/Tara_metatranscriptomic"
reference = list(csv.reader(open(f"{base1}/OM-RGC_v2.tsv"), delimiter='\t'))

def get_id_v2(ID_list):
	gene_id = set()
	with open(ID_list, 'r') as f:
		all_id = f.readlines()
	print(all_id[:10])
	for id in all_id:
		id = id.strip()
		gene_id.add(id)
	return gene_id

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
	# with open(f"{gene_name}_OMRGC_ID_length.csv", 'w') as f:f.write("\n".join(out))
	return gene_length_dict

def get_gene_specific_ref_through_COG(search_string):
	cog_out = set()
	cogfile = open(f'../COG_numbers_categories_annotations.txt', 'r')
	for line in cogfile:
		columns = line.split('\t')
		cog_number, cog_category = columns[0], columns[2]
		if search_string in cog_category:
			cog_out.add(cog_number)
	cogfile.close()
	out = []
	gene_length_dict = {}
	for line in reference:
		if line[3] in cog_out:
			gene_id, COG, gene_length = line[1], line[3], len(line[13])
			out.append(f"{gene_id}, {COG}, {gene_length}")
			gene_length_dict[gene_id] = gene_length
	# with open(f"{search_string}_OMRGC_ID_length.csv", 'w') as f:f.write("\n".join(out))
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
	with open(f"{base1}/{gene_name}_{DNA_or_RNA}_coverage_subset.csv", 'w') as f:f.write("\n".join(out))

def cal_sum_num_reads_mapped(gene_name, DNA_or_RNA):
	counter = 0
	subset = list(csv.reader(open(f"{base1}/{gene_name}_{DNA_or_RNA}_coverage_subset.csv") , delimiter='\t'))
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
	normalization_factor = pd.read_csv(f"{base1}/sample_readCounts_readLength_{DNA_or_RNA}.tsv", sep = "\t")
	sum_len = sum_len.merge(normalization_factor, on="sample", how="inner")
	print(sum_len)
	sum_len['read_counts'] = sum_len['sum_bases']/sum_len['Avg_read_len']
	sum_len['prop_reads_mapped_back_genes'] = sum_len['read_counts']/sum_len['HQ_reads']
	sum_len.to_csv(f'{base1}/prop_{DNA_or_RNA}_reads_mapped_back_{gene_name}.csv', sep=',', index=False)
	return sum_len

def serialize(gene_length_dict, gene_name):
	subset_coverage_rows(gene_length_dict, gene_name, f"{base1}/OM-RGC_v2_gene_profile_metaG.tsv", "DNA")
	DNA_category_sum = cal_sum_num_reads_mapped(gene_name, "DNA")
	DNA_category_prop = cal_prop_reads_mapped(DNA_category_sum, gene_name, "DNA")
	print(DNA_category_prop)
	subset_coverage_rows(gene_length_dict, gene_name, f"{base2}/OM-RGC_v2_gene_profile_metaT.tsv", "RNA")
	RNA_category_sum = cal_sum_num_reads_mapped(gene_name, "RNA")
	RNA_category_prop = cal_prop_reads_mapped(RNA_category_sum, gene_name, "RNA")
	print(RNA_category_prop)

def serialized_COG_category(COG_category):
	COG_length_dict = get_gene_specific_ref_through_COG(COG_category)
	print(f"done with dict, there are {len(COG_length_dict)} COG genes")
	serialize(COG_length_dict, COG_category)

def serialized_genes(gene_id_set, gene_name):
	gene_length_dict = get_gene_specific_ref(gene_id_set, gene_name)
	print(f"done with dict, there are {len(gene_length_dict)} genes")
	serialize(gene_length_dict, gene_name)

def merge_abun_helper(gene, DNA_or_RNA):
	base = pd.read_csv(f"{base1}/prop_{DNA_or_RNA}_reads_mapped_back_{gene}.csv")
	base = base[["sample", "prop_reads_mapped_back_genes"]]
	gene = gene.replace("transport_and_metabolism", "TM")
	base.columns = [f"connector_{DNA_or_RNA}", f"{DNA_or_RNA}_{gene}"]
	return base

def merge_abundance(all_genes, DNA_or_RNA):
	base = merge_abun_helper(all_genes[0], DNA_or_RNA)
	for gene in all_genes[1:]:
		appending = merge_abun_helper(gene, DNA_or_RNA)
		base = base.merge(appending, on = f"connector_{DNA_or_RNA}")
	return base

# biofilm_id_set = get_id("OM-RGC_v2_Reference_Pieces/Biofilm_BLAST_Result/merge_Biofilm_hits_no_BapA.txt")
# biofilm_length_dict = get_gene_specific_ref(biofilm_id_set, "biofilm")

trans_id_set = get_id(f"{base1}/OM-RGC_v2_Reference_Pieces/Transposase_BLAST_Result/merged_Transposase_hits.txt")
CAZ_id_set = get_id_v2(f"{base1}/all_CAZyme_norm.txt")
pep_id_set = get_id_v2(f"{base1}/all_peptidase_norm.txt")
sect_CAZ_id_set = get_id_v2(f"{base1}/secretory_CAZyme.txt")
sect_pep_id_set = get_id_v2(f"{base1}/secretory_peptidase.txt")

single_genes = ["transposase", "CAZyme", "peptidase", "secretory_CAZyme", "secretory_peptidase"]
gene_ide_set = [trans_id_set, CAZ_id_set, pep_id_set, sect_CAZ_id_set, sect_pep_id_set]

for i in range(len(single_genes)):
	print("doing: " + single_genes[i])
	serialized_genes(gene_ide_set[i], single_genes[i])

serialized_genes(CAZ_id_set, "CAZyme")
serialized_genes(pep_id_set, "peptidase")

COG_list = ["Lipid_transport_and_metabolism", "Coenzyme_transport_and_metabolism", "Signal_transduction_mechanisms", 
"Defense_mechanisms", "Energy_production_and_conversion", "Replication_recombination_and_repair", "Cell_wall",
"Amino_acid_transport_and_metabolism","Cell_motility", "Extracellular_struct", "Cytoskeleton",
"Carbohydrate_transport_and_metabolism", "Inorganic_ion_transport_and_metabolism", "Nucleotide_transport_and_metabolism"]

for COG_category in COG_list:
	print("doing: " + COG_category)
	serialized_COG_category(COG_category)

all_genes = single_genes + COG_list
DNA_merged = merge_abundance(all_genes, "DNA")
DNA_merged.to_csv(f"../Figure Generating/data/tara_DNA_genes_abundance.csv", index=False)
RNA_merged = merge_abundance(all_genes, "RNA")
RNA_merged.to_csv(f"../Figure Generating/data/tara_RNA_genes_abundance.csv", index=False)

# DNA
# subset_coverage_rows(biofilm_length_dict, "biofilm", "OM-RGC_v2_gene_profile_metaG.tsv", "DNA")
# subset_coverage_rows(trans_length_dict, "transposase", "OM-RGC_v2_gene_profile_metaG.tsv", "DNA")
# biofilm_sum = cal_sum_num_reads_mapped("biofilm", "DNA")
# biofilm_prop = cal_prop_reads_mapped(biofilm_sum, "biofilm", "DNA")
# transposase_sum = cal_sum_num_reads_mapped("transposase", "DNA")
# transposase_prop = cal_prop_reads_mapped(transposase_sum, "transposase", "DNA")

# RNA
# subset_coverage_rows(biofilm_length_dict, "biofilm", "../Tara_metatranscriptomic/OM-RGC_v2_gene_profile_metaT.tsv", "RNA")
# subset_coverage_rows(trans_length_dict, "transposase", "../Tara_metatranscriptomic/OM-RGC_v2_gene_profile_metaT.tsv", "RNA")
# biofilm_sum = cal_sum_num_reads_mapped("biofilm", "RNA")
# biofilm_prop = cal_prop_reads_mapped(biofilm_sum, "biofilm", "RNA")
# transposase_sum = cal_sum_num_reads_mapped("transposase", "RNA")
# transposase_prop = cal_prop_reads_mapped(transposase_sum, "transposase", "RNA")

