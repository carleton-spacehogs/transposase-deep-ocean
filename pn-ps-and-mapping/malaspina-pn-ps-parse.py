# Written by Jimmy Zhong (zhongj2@carleton.edu), Carleton '23 under Professor Rika Anderson; date: July 16 2021
import glob
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import csv
import time
import random

integrons = list(csv.reader(open("malaspina-deep-co-assembly.integrons"), delimiter='\t'))[1:] 
pn_ps = list(csv.reader(open("pNpS.txt"), delimiter='\t'))[1:]
transposases = list(csv.reader(open("cog_function_transposase.txt"), delimiter='\t'))[1:]
biofilms = list(csv.reader(open("cog_function_biofilm.txt"), delimiter='\t'))[1:]
all_gene_calls = list(csv.reader(open("EVERYTHING-gene_calls.txt"), delimiter='\t'))[1:]
# prokka_file = list(csv.reader(open("../malaspina-deep-co-assembly_prokka/deep_cleaned.gff"), delimiter='\t'))

print("file reading done!")

def within100(start1, start2, end1, end2):
	if abs(int(start1) - int(start2)) + abs(int(end1) - int(end2)) < 200:
		return True
	else:
		return False

def get_gene_category(gene_caller_id):
	if gene_caller_id in interest_gene_ids:
		return interest_gene_ids[gene_caller_id]
	else:
		return "Normal"

# def parse_prokka_line(prokka_ID):
# 	this_line_dict = {}
# 	prokka_ID = prokka_ID.strip("\n")
# 	categories = prokka_ID.split(";")
# 	for each in categories:
# 		if "=" in each:
# 			key, value = each.split("=")
# 			this_line_dict[key] = value
# 	return this_line_dict

# def get_COG_functional_group(gene_info_dict, COG_dict):
# 	COG_fun_group = 'NA'
# 	COG_num = -1
# 	if "db_xref" in gene_info_dict:
# 		COG_num = gene_info_dict["db_xref"].replace("COG:","")
# 		if COG_num in COG_dict:
# 			COG_fun_group = COG_dict[COG_num]
# 	return str(COG_num), str(COG_fun_group)

# def COG_dict():
# 	#first, make a dictionary matching the COG number to the COG category.
# 	D = {}
# 	cogfile = open('/usr/local/data/Python_scripts/COG_numbers_categories_annotations.txt', 'r')
# 	for line in cogfile:
# 		columns = line.split('\t')
# 		cog_number = columns[0]
# 		cog_category = columns[2].replace("_and", "")
# 		# some error fixing
# 		if cog_category == "Amino_acid_transport_and_metabolism/Nucletide_transport_metabolism": 
# 			cog_category = "Amino_acid_transport_and_metabolism/Nucleotide_transport_metabolism"
# 		if cog_number in D: # another function group!
# 			D[cog_number] = D[cog_number] + " " + columns[0]
# 		else:
# 			D[cog_number] = cog_category
# 	cogfile.close()
# 	return D

interestings = [transposases, biofilms]
interest_text = ["standalone_transposase", "standalone_biofilm"]
interest_gene_ids = {}
for i in range(len(interestings)):
	name = interest_text[i]
	for line in interestings[i]:
		gene_caller_id = line[0]
		interest_gene_ids[gene_caller_id] = name

pn_ps_cur = 0
out = []
for i in range(len(all_gene_calls)):
	gene_call_id = all_gene_calls[i][0]
	if gene_call_id == pn_ps[pn_ps_cur][0]:
		pn_ps_ratio = pn_ps[pn_ps_cur][3]
		gene_category = get_gene_category(gene_call_id)
		out_line = [gene_call_id, pn_ps_ratio, gene_category] + all_gene_calls[i][1:-1]
		out.append(out_line)
		pn_ps_cur += 1
	if pn_ps_cur%200000 == 0: print("I am at line :" + str(pn_ps_cur))
	if pn_ps_cur == len(pn_ps): break
	if not i<pn_ps_cur: print("ERROR!!! some caller id in pNpS.txt is not in EVERYTHING-gene_calls.txt")

# out_cur = 0
# COG_dict = COG_dict()
# for i in range(len(prokka_file)):
# 	if prokka_file[i][0] == out[out_cur][3]: # contig name matches
# 		prokka_start, prokka_end = int(prokka_file[i][3]), int(prokka_file[i][4])
# 		anvi_start, anvi_end = int(out[out_cur][4]), int(out[out_cur][5])
# 		if within100(anvi_start, prokka_start, anvi_end, prokka_end):
# 			gene_info_dict = parse_prokka_line(prokka_file[i][8]) # exmaple: ID=DMEEPCGM_00001;Name=cusA_1;dbxref=COG:COG3696;gene=cusA_1;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:P38054;locus_tag=DMEEPCGM_00001;product=Cation efflux system protein CusA
# 			out[out_cur][13] = gene_info_dict["product"] # product == gene function
# 			out[out_cur][14], out[out_cur][15] = get_COG_functional_group(gene_info_dict, COG_dict)# COG_num, COG_functional_group
# 			out_cur += 1
# 	if out_cur%200000 == 0: print("I am at line :" + str(pn_ps_cur))
# 	if out_cur == len(out): break
# 	if not i>out_cur: print("ERROR!!! in the second step")

print("done merging")

integron_contigs = set()
for line in integrons:
	integron_contigs.add(line[1])

found = 0
for i in range(len(out)):
	contig_name = out[i][3]
	if contig_name in integron_contigs:
		for line in integrons:
			if contig_name == line[1]:
				match = 0
				if within100(out[i][4], line[3], out[i][5], line[4]): # anvi and integron start end respetively
					match = True
					found += 1
					gene_category = line[8] + "_" + line[10]
					if "toxin" in out[i][11]: gene_category = line[10] + "_toxin"
					if out[i][2] == "standalone_transposase": gene_category = line[10] + "_transposase"
					if out[i][2] == "standalone_biofilm": gene_category = line[10] + "_biofilm"
					out[i][2] = gene_category
				if not match:
					print(str(out[i]) + " integron_finder calls: " + str(line[3:5]))
print("I found this many integrons: " + str(found))


# integrons = []
# for line in out:
# 	if not(line[2] == "Normal" or line[2] == "standalone_transposase" or line[2] == "standalone_biofilm"):
# 		integrons.append(line)


# with open("pNpS_merged_subsampled.csv", 'w') as f:
# 	write = csv.writer(f)
# 	write.writerow(out_header)
# 	write.writerows(subsampled)
# print("subsample writing file done")

header = ["gene_caller_id", "pn_ps", "gene_category", "contig", "start", "stop", "direction", 
"COG20_CATEGORY", "COG20_CATEGORY (ACCESSION)", "COG20_PATHWAY", "COG20_PATHWAY (ACCESSION)", 
"COG20_FUNCTION", "COG20_FUNCTION (ACCESSION)"]

# with open("pNpS_integron.csv", 'w') as f:
# 	write = csv.writer(f)
# 	write.writerow(header)
# 	write.writerows(integrons)
# print("integron writing file done")

with open("pNpS_merged.csv", 'w') as f:
	write = csv.writer(f)
	write.writerow(header)
	write.writerows(out)

print("full file writing done")

df = pd.read_csv("pNpS_merged.csv")
df = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)]
print(df[["gene_category", 'pn_ps']].groupby("gene_category").describe())
print(df[["gene_category", 'pn_ps']].groupby("gene_category").mean(numeric_only=True))
print(df[["gene_category", 'pn_ps']].groupby("gene_category").median(numeric_only=True))

