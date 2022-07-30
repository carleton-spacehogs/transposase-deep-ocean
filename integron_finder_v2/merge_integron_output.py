#!/bin/python3
from email.charset import add_charset
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
deep_rt = "/researchdrive/zhongj2/deep_ocean_bins/"

def get_anvi_db(root, ocean):
	gene_call_f = f"{root}/{ocean}_bins_v2/all-gene-calls.txt"
	if ocean == "deep":
		gene_call_f = f"{deep_rt}/anvio-gene-call.txt"
	all_gene_calls = pd.read_csv(gene_call_f, sep='\t')
	all_gene_calls.drop(columns = ['source','version'], inplace = True)
	return all_gene_calls

def get_anvi_COG(root, ocean):
	COG_file = f"{root}/{ocean}_bins_v2/all-COG-categories.txt"
	if ocean == "deep":
		COG_file = f"{deep_rt}/all-anvi-cog-category.tsv"
	category = pd.read_csv(COG_file, sep="\t")
	category.drop(columns = ['source','e_value'], inplace = True)
	return category

o_depth = ocean_depths()
oceans = list(o_depth.keys()) + ["deep"]
integron_dir="/researchdrive/zhongj2/integron_finder_tara_contig"
anvio_db_dir="/researchdrive/zhongj2/Tara_Oceans_Bins"

merged_integrons = pd.DataFrame()
for ocean in oceans:
	print(ocean)
	integron_dict, int_f = {}, f"{integron_dir}/integron-V2-{ocean}.txt"
	integron_dict = integron_dict_helper(int_f, integron_dict)
	anvi_contig_gene_caller_df = get_anvi_db(anvio_db_dir, ocean)
	anvi_COG_category = get_anvi_COG(anvio_db_dir, ocean)
	integrons = merge_anvi_gene_call(anvi_contig_gene_caller_df, ocean, integron_dict)
	integrons = integrons.merge(anvi_COG_category, on="gene_callers_id", how = "left")
	integrons = explode_COG(integrons, 10)
	merged_integrons = merged_integrons.append(integrons)

merged_integrons = merged_integrons[merged_integrons.integron.str.contains("CALIN|complete")]

functional_integrons = merged_integrons[merged_integrons.function != "nan"]
functional_integrons.to_csv("metagenome-integron-with-COG-functions.csv",sep=',', index=False)


# pattern = "nan|Function unknown|General function prediction only"
# merged_integrons[~merged_integrons.function.str.contains(pattern)] # 7966
# merged_integrons[merged_integrons.function.str.contains("Defense mechanisms")] # 734
