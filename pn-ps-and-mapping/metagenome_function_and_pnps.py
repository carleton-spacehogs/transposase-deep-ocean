# Written by Jimmy Zhong (zhongj2@carleton.edu), Carleton '23 under Professor Rika Anderson; date: August 31nd 2021
'''
cd /workspace/data/zhongj/Transposase_Project/integron_finder_tara_contig
python3
exec(open('find_cassette_function_and_pnps_august24.py').read())
'''
import glob
import os
import subprocess
import sys
import csv
import time
import random
import pandas as pd
import numpy as np
from collections import Counter
ocean_list = ["ARS", "CPC", "EAC", "IN", "NAT", "MED", "NP", "RS", "SAT", "SP"]

# mv tara_arabian_SECONDARY_contigs.integrons ARS.integrons
# mv tara_chileperucoastal_SECONDARY_contigs.integrons CPC.integrons
# mv tara_eafricacoastal_SECONDARY_contigs.integrons EAC.integrons
# mv tara_indianmonsoon_SECONDARY_contigs.integrons IN.integrons
# mv tara_mediterranean_SECONDARY_contigs.integrons MED.integrons
# mv tara_northatlantic_SECONDARY_contigs.integrons NAT.integrons
# mv tara_northpacific_SECONDARY_contigs.integrons NP.integrons
# mv tara_redsea_SECONDARY_contigs.integrons RS.integrons
# mv tara_southatlantic_SECONDARY_contigs.integrons SAT.integrons
# mv tara_southpacific_SECONDARY_contigs.integrons SP.integrons

integron_merged_colnames = ['gene_callers_id', 'contig', 'start', 'stop', 'direction', 'partial', 'call_type',
'pNpS', 'accession', 'gene_function', 'e_value', 'category', 'element', 'cassette_begin', 'cassette_end',
'type_molecule', 'from_integron_type']

def within(start1, start2, end1, end2):
	if abs(int(start1)-int(start2)) < 100 and abs(int(end1)-int(end2)) < 100:
		return True
	else:
		return False

def parse_anvio_files(ocean):
	all_gene_calls = pd.read_csv(f"../Tara_Oceans_Bins/{ocean}_bins_v2/all-gene-calls.txt", sep='\t')
	all_gene_calls.drop(columns = ['source','version'], inplace = True)
	func = pd.read_csv(f"../Tara_Oceans_Bins/{ocean}_bins_v2/all-COG-functions.txt", sep='\t')
	func.drop(columns = ['source'], inplace = True)
	category = pd.read_csv(f"../Tara_Oceans_Bins/{ocean}_bins_v2/all-COG-categories.txt", sep="\t")
	category.drop(columns = ['source','e_value', 'accession'], inplace = True)
	pnps_file = f"../Tara_Oceans_Bins/{ocean}_bins_v2/DNA-pn-ps/pNpS.txt" # I only compute integron contigs for EAC, SP, and SAT. Computational restriction, not enough memory
	if not os.path.isfile(pnps_file): 
		pnps_file = f"../Tara_Oceans_Bins/{ocean}_bins_v2/integron-pn-ps/pNpS.txt"
	pnps = pd.read_csv(pnps_file, sep='\t')
	pnps.drop(columns = ["Unnamed: 0", "sample_id"], inplace = True)
	pnps.rename(columns={"corresponding_gene_call":"gene_callers_id", "pNpS_gene_reference":"pNpS"}, inplace=True)
	print(pnps)
	print(all_gene_calls)
	all_gene_calls = all_gene_calls.merge(pnps, on="gene_callers_id", how='left')
	gene_call_func_category = (all_gene_calls.merge(func, how='outer', on="gene_callers_id")).merge(category, on="gene_callers_id", how='outer')
	gene_call_func_category.rename(columns={"function_y": "category", "function_x": "gene_function"}, inplace = True)
	return gene_call_func_category.values.tolist()

def parse_integrons(ocean):
	with open(f"{ocean}.integrons") as file: 
		hits = file.readlines()[2:]
	integrons = []
	# ["contig,element,pos_begin,pos_end,type_molecule,from_integron_type\n"]
	# all anvio's contig names are c_1, c_2, etc
	for line in hits:
		line = line.split("\t")
		if len(line) > 10: integrons.append(["c_" + line[1], line[2], line[3], line[4], line[7], line[10]])
	return integrons

for ocean in ocean_list:
	integrons = parse_integrons(ocean)
	integrons_set = set(integron[0] for integron in integrons)
	gene_call_func_category_list = parse_anvio_files(ocean)
	out = []
	for anvio_gene_call in gene_call_func_category_list:
		contig_name = anvio_gene_call[1]
		if contig_name in integrons_set: # check whether it's one of the integron containing contigs
			gene_call_start, gene_call_end = anvio_gene_call[2], anvio_gene_call[3]
			for cassette in integrons:
				if contig_name == cassette[0]:
					cassette_start, cassette_end = cassette[2], cassette[3]
					if within(cassette_start, gene_call_start, cassette_end, gene_call_end): 
						out.append(anvio_gene_call + cassette[1:])
	merged_integrons = pd.DataFrame(out, columns = integron_merged_colnames)
	merged_integrons.to_csv(path_or_buf=f'{ocean}_merged_integrons.csv', sep=',', index=False)

def cat_all_ocean_integron(ocean_list):
	catted = pd.DataFrame(columns = integron_merged_colnames)
	for ocean in ocean_list:
		additional = pd.read_csv(f"{ocean}_merged_integrons.csv")
		additional['ocean']=ocean
		print(ocean + ", length: " + str(len(additional.index)))
		catted = catted.append(additional)
	return catted

def cat_all_ocean_random(ocean_list):
	catted = pd.DataFrame()
	for ocean in ocean_list:
		normal_all = pd.read_csv(f"../Tara_Oceans_Bins/{ocean}_bins_v2/all-COG-categories.txt", sep='\t').rename(columns={"function": "category"})
		catted = catted.append(normal_all)
	return catted

normal_gene_call = cat_all_ocean_random(ocean_list)

all_ocean_integron = cat_all_ocean_integron(ocean_list)
all_ocean_integron.to_csv(path_or_buf='all_ocean_merged_integrons.csv', sep=',', index=False)

# filter out the integrase of the cassette seequence, ~ negates condition
all_ocean_integron = all_ocean_integron[~all_ocean_integron['from_integron_type'].str.contains("In")]
all_ocean_integron = all_ocean_integron[~all_ocean_integron['type_molecule'].str.contains("attC")]

def count_cog_proprotion(pd_dataframe):
	tmp = pd_dataframe[pd_dataframe.category.notnull()]
	# if 'from_integron_type' in tmp: 
	# 	tmp = tmp[~tmp['from_integron_type'].str.contains("In")] 
	func_count_list = tmp['category'].tolist()
	func_count_list = [word for line in func_count_list for word in line.split("!!!")] # split !!! into 2 different functions
	func_Counter = Counter(func_count_list)
	out, total = [], len(func_count_list)
	for tag, count in func_Counter.items():  
		out.append([str(tag), str(count), str(count/total)])
	return out, total

integron_count, total_integron = count_cog_proprotion(all_ocean_integron)
# SAT_gene_call = pd.read_csv(f"../Tara_Oceans_Bins/SAT_bins_v2/all-COG-categories.txt", sep='\t').rename(columns={"function": "category"})
integron_count, total_integron = count_cog_proprotion(all_ocean_integron)
normal_count, total_normal = count_cog_proprotion(normal_gene_call)

integron_df = pd.DataFrame(integron_count, columns = ['COG_function','integron_count','integron_prop'])
normal_df = pd.DataFrame(normal_count, columns = ['COG_function','normal_count','normal_prop'])
merged = integron_df.merge(normal_df, on="COG_function", how='inner')
merged.sort_values(by = "COG_function", axis=0)
merged.loc[len(merged.index)] = ['~Total', total_integron, "1", total_normal, "1"] 

merged.to_csv(path_or_buf="all_integron_func_category_count.csv", sep=',', index=False)
 
all_integron = pd.read_csv("all_ocean_merged_integrons.csv")
all_integron = all_integron[["pNpS", "category", 'type_molecule', 'from_integron_type']]
non_null_integron = all_integron.dropna()
non_null_integron.to_csv(path_or_buf="integron_pnps_clean.csv", sep=',', index=False)

# pd.options.display.max_colwidth = 100
# pd.options.display.max_rows = 50

# non_null_integron['category']=non_null_integron['category'].str.split('!!!')
# non_null_integron.explode('category')

print("finish the first half, now calculate the mean bin pnps for all bins")

pnps_bin_contig = pd.DataFrame(columns = ["gene_callers_id","contig","pNpS","bin_name"])
for ocean in ocean_list:
	all_gene_calls = pd.read_csv(f"./{ocean}_bins_v2/all-gene-calls.txt", sep='\t')
	all_gene_calls = all_gene_calls[["gene_callers_id", "contig"]]
	bin_pnps_file = f"./{ocean}_bins_v2/DNA-pn-ps/pNpS.txt" # I only compute integron contigs for EAC, SP, and SAT. Computational restriction, not enough memory
	if not os.path.isfile(bin_pnps_file): bin_pnps_file = bin_pnps_file.replace("DNA", "bin")
	pnps = pd.read_csv(bin_pnps_file, sep='\t') # pnps=pd.read_csv("bin-pn-ps/pNpS.txt", sep='\t')
	pnps.drop(columns = ["Unnamed: 0", "sample_id"], inplace = True)
	pnps.rename(columns={"corresponding_gene_call":"gene_callers_id", "pNpS_gene_reference":"pNpS"}, inplace=True)
	tmp = all_gene_calls.merge(pnps, how="inner", on="gene_callers_id")
	binning_info = pd.read_csv(f"./{ocean}_bins_v2/{ocean}_binning_info_anvio.txt", names=["contig", "bin_name"], sep='\t')
	tmp = tmp.merge(binning_info, on="contig")
	func = pd.read_csv(f"./{ocean}_bins_v2/all-COG-functions.txt", sep='\t')
	func.drop(columns = ['source'], inplace = True)
	tmp = tmp.merge(func, how='outer', on="gene_callers_id")
	tmp = tmp[~tmp.isin([np.nan, np.inf, -np.inf]).any(1)] # filter out inf, -inf, and nan
	pnps_bin_contig = pnps_bin_contig.append(tmp)

# # get bins' median pnps
# bin_median_pnps = pnps_bin_contig.groupby(["bin_name"])['pNpS'].agg(mean_bin_pnps='mean', median_bin_pnps='median', count='size').reset_index()

# df_mask=bin_median_pnps['count']>=10
# bin_median_pnps = bin_median_pnps[df_mask]

# bin_median_pnps.rename(columns={"bin_name":"bin"}, inplace=True)
# bin_median_pnps.to_csv('bin_median_pnps.csv', sep=',', index=False)


# show that bins with high pnps does not also have a ton of transposases
all_integrons = pd.read_csv('all_ocean_merged_integrons.csv')

pnps_bin_contig = pd.DataFrame(columns = ["gene_callers_id","contig","pNpS","bin_name"])
for ocean in ocean_list:
	integrons = parse_integrons(ocean)
	integrons_set = set(integron[0] for integron in integrons)
	all_gene_calls = pd.read_csv(f"./{ocean}_bins_v2/all-gene-calls.txt", sep='\t')
	all_gene_calls = all_gene_calls[["gene_callers_id", "contig"]]
	bin_pnps_file = f"./{ocean}_bins_v2/DNA-pn-ps/pNpS.txt" # I only compute integron contigs for EAC, SP, and SAT. Computational restriction, not enough memory
	if not os.path.isfile(bin_pnps_file): bin_pnps_file = bin_pnps_file.replace("DNA", "bin")
	pnps = pd.read_csv(bin_pnps_file, sep='\t') # pnps=pd.read_csv("bin-pn-ps/pNpS.txt", sep='\t')
	pnps.drop(columns = ["Unnamed: 0", "sample_id"], inplace = True)
	pnps.rename(columns={"corresponding_gene_call":"gene_callers_id", "pNpS_gene_reference":"pNpS"}, inplace=True)
	tmp = all_gene_calls.merge(pnps, how="inner", on="gene_callers_id")
	binning_info = pd.read_csv(f"./{ocean}_bins_v2/{ocean}_binning_info_anvio.txt", names=["contig", "bin_name"], sep='\t')
	tmp = tmp.merge(binning_info, on="contig")
	tmp = tmp[~tmp.isin([np.nan, np.inf, -np.inf]).any(1)] # filter out inf, -inf, and nan
	pnps_bin_contig = pnps_bin_contig.append(tmp)

bin_median_pnps = pnps_bin_contig.groupby(["bin_name"])['pNpS'].agg(mean_bin_pnps='mean', median_bin_pnps='median', count='size').reset_index()

df_mask=bin_median_pnps['count']>=10
bin_median_pnps = bin_median_pnps[df_mask]

bin_median_pnps.rename(columns={"bin_name":"bin"}, inplace=True)
bin_median_pnps.to_csv('bin_median_pnps.csv', sep=',', index=False)

def get_upper_range(df):
	Q3, Q1 = np.quantile(df['pNpS'], 0.75), np.quantile(df['pNpS'], 0.25)
	IQR = Q3 - Q1
	return 1.5*IQR + Q3
upper_range = get_upper_range(pnps_bin_contig)

# get pnps for all transposases in bins
transposase_in_bins = pnps_bin_contig[pnps_bin_contig['function'].str.contains("transposase")]
non_transposase = pnps_bin_contig[~pnps_bin_contig['function'].str.contains("transposase")]
print(transposase_in_bins.query('pNpS < @upper_range').pNpS.describe())

print(scipy.stats.ttest_ind(non_transposase.query('pNpS < @upper_range').pNpS, transposase_in_bins.query('pNpS < @upper_range').pNpS))
