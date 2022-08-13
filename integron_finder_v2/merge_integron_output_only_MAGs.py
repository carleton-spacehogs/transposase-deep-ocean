#!/bin/python3
import sys
import pandas as pd
sys.path.insert(0, '../MAG_pnps_redux')
from MAG_integron_selection import explode_COG, integron_dict_helper, merge_anvi_gene_call
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths

''' Merge the integron finder output from all oceans, and assign the COG functional
category for them.

Written by Jimmy Zhong at Carleton College, MN, working under Professor Rika Anderson's
Spacehog Bioinformaitcs lab, July 30th 2022.
'''
MAG_dir = "/researchdrive/zhongj2/MAG_pnps_redux"
o_depth = ocean_depths()
oceans = list(o_depth.keys()) + ["deep"]

def get_anvi_db(root, ocean):
	gene_call_f = f"{MAG_dir}/{ocean}/all_bins_db/gene_callers_id-contig_full.txt"
	all_gene_calls = pd.read_csv(gene_call_f, sep='\t')
	all_gene_calls.drop(columns = ['source','version'], inplace = True)
	return all_gene_calls

def get_anvi_COG(root, ocean):
	COG_file = f"{root}/{ocean}/all_bins_db/anvi_genes_COG_categories.tsv"
	category = pd.read_csv(COG_file, sep="\t")
	category.drop(columns = ['e_value'], inplace = True)
	return category

merged_integrons = pd.DataFrame()
for ocean in oceans:
	print(ocean)
	int_f = f"{MAG_dir}/{ocean}/all_bins_db/Results_Integron_Finder_{ocean}_all_bins/all_integron.txt"
	integron_dict = {}
	integron_dict = integron_dict_helper(int_f, integron_dict)
	anvi_contig_gene_caller_df = get_anvi_db(MAG_dir, ocean)
	anvi_COG_category = get_anvi_COG(MAG_dir, ocean)
	integrons = merge_anvi_gene_call(anvi_contig_gene_caller_df, ocean, integron_dict)
	integrons = integrons.merge(anvi_COG_category, on="gene_callers_id", how = "left")
	integrons = explode_COG(integrons, 10)
	merged_integrons = merged_integrons.append(integrons)

merged_integrons = merged_integrons[merged_integrons.integron.str.contains("CALIN|complete")]
merged_integrons.to_csv("MAG-integron-all.csv",sep=',', index=False)

functional_integrons = merged_integrons[merged_integrons.function != "nan"]
functional_integrons.to_csv("metagenome-integron-with-COG-functions.csv",sep=',', index=False)

