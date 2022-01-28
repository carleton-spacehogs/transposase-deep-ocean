# Written by Jimmy Zhong (zhongj2@carleton.edu), Carleton '23 under Professor Rika Anderson; date: Jan 27 2022
import glob
import os
import sys
import csv
import pandas as pd

def dict_to_list(dict):
	rows = []
	for key, value in dict.items(): rows.append([key, value])
	return rows

def target_id_set(pNpS2csv):
	target_id = set()
	ids = pd.read_csv(pNpS2csv, usecols=["gene_caller_id"]).values.tolist()
	for i in ids: target_id.add(i[0])
	return target_id

integron_caller_id = target_id_set("./pNpS2_integron.csv")
normal_caller_id = target_id_set("./pNpS2_non_integrons_subsampled.csv")
pd_anvio = pd.read_csv(f"./EVERYTHING-gene_calls.txt", sep='\t', usecols=["gene_callers_id", "COG20_CATEGORY"])
pd_anvio = pd_anvio[pd_anvio.COG20_CATEGORY.notnull()]
all_gene_calls = pd_anvio.values.tolist()

print("done reading files")

def get_function_group(in_dict):
	func_dict = {}
	for line in all_gene_calls:
		gene_call_id, func = line[0], line[1]
		if gene_call_id in in_dict:
			funcs = func.split("!!!")
			for func in funcs:
				if func in func_dict:
					func_dict[func] += 1
				else:
					func_dict[func] = 1
	# func_dict_group = dict_to_list(func_dict)
	return func_dict

print("get integron functional groups")
integron_group = get_function_group(integron_caller_id)

print("now getting sampled normal genes")
normal_group = get_function_group(normal_caller_id)

integron_total, normal_total = 0, 0 # some functional groups are not present in integron_group
funcs = list(integron_group.keys())
funcs.sort()
for func in funcs:
	normal_total += normal_group[func]
	integron_total += integron_group[func]

merge_out = []
header = ["COG_func_group", "count_all", "count_integron", "proportion_all", "proportion_integron"]
for func in funcs:
	nor_c, int_c = normal_group[func], integron_group[func] # normal, integron count of this funcion
	to_add = [func, nor_c, int_c, nor_c/normal_total, int_c/integron_total]
	merge_out.append(to_add)

with open("merged_function_group_proportion.csv", 'w') as f:
	write = csv.writer(f)
	write.writerow(header)
	write.writerows(merge_out)
