# Written by Jimmy Zhong (zhongj2@carleton.edu), Carleton '23 under Professor Rika Anderson; date: Jan 27 2022
import glob
import pandas as pd
import numpy as np
import os
import subprocess
import sys
import csv
import time
import random
import math

def parse_anvio_files():
	all_gene_calls = pd.read_csv(f"./EVERYTHING-gene_calls.txt", sep='\t',
	usecols=["gene_callers_id", "contig", "start", "stop", "COG20_CATEGORY (ACCESSION)"])
	pnps = pd.read_csv("pNpS.txt", sep='\t')
	pnps.drop(columns = ["Unnamed: 0", "sample_id"], inplace = True)
	pnps.rename(columns={"corresponding_gene_call":"gene_callers_id", "pNpS_gene_reference":"pnps"}, inplace=True)
	all_gene_calls = all_gene_calls.merge(pnps, on="gene_callers_id", how='left')
	return all_gene_calls

def target_id_set(blastp_file):
	target_id = set()
	with open(blastp_file, newline='') as f:
		blast_result=list(csv.reader(f,delimiter='\t'))
	for line in blast_result:
		target_id.add(line[1])
	return target_id

def within100(start1, start2, end1, end2):
	if abs(int(start1) - int(start2)) + abs(int(end1) - int(end2)) < 100:
		return True
	else:
		return False

pd_anvio = parse_anvio_files()
gene_calls = pd_anvio.values.tolist()
transposase = target_id_set("../deep_metagenome_transposase_BLAST/deep_anvi-genecall_transposase.blastp")
biofilm = target_id_set("../deep_metagenome_transposase_BLAST/deep_anvi-genecall_biofilm.blastp")
integrons = pd.read_csv("malaspina-deep-co-assembly.integrons", delimiter='\t')

print("file reading done!")

integrons.drop(columns = ['ID_integron','element','strand','evalue','model',
'annotation','default','distance_2attC',"considered_topology"], inplace = True)
cassette = integrons[integrons['type_elt']=="protein"].values.tolist()
cassette_contigs = set()
for line in cassette: cassette_contigs.add(line[0])

def isNaN(num):
	return num!= num

header = ["gene_caller_id", "COG20_CATEGORY", "pnps", "gene_type"]
integrons = []
non_integrons = []
for gene_call in gene_calls:
	# usecols=["gene_callers_id", "contig", "start", "stop", "COG20_CATEGORY (ACCESSION)"]) + pnps
	contig_name = gene_call[1]
	gene_caller_id, func, pnps = gene_call[0], gene_call[-2], gene_call[-1]
	new_row = [gene_caller_id, func, pnps, "__"] 
	if contig_name in cassette_contigs:
		for line in cassette:
			if contig_name == line[0]:
				gc_s, gc_e = gene_call[2], gene_call[3] #gene call start, end
				c_s, c_e = line[1], line[2] #cassette start, end
				if within100(gc_s, c_s, gc_e, c_e):
					if isNaN(func): new_row[-1] = "C_" # cassette no gene call
					elif "V" in func: new_row[-1] = "CD" # defense mechanisms
					else: new_row[-1] = "CN" # non-defense gene calls
					integrons.append(new_row)
	else:
		if not isNaN(pnps): # do not keep non-cassette rows without pnps
			if gene_caller_id in transposase:
				new_row[-1] = "_T" # transposase
			elif gene_caller_id in biofilm:
				new_row[-1] = "_B" # biofilm
			elif random.randint(0, 7) == 0:
				non_integrons.append(new_row)

print("I found this many integrons: ", len(integrons))
print("I found this many random: ", len(non_integrons))


with open("pNpS2_non_integrons_subsampled.csv", 'w') as f:
	write = csv.writer(f)
	write.writerow(header)
	write.writerows(non_integrons)

with open("pNpS2_integron.csv", 'w') as f:
	write = csv.writer(f)
	write.writerow(header)
	write.writerows(integrons)

# with open("pNpS_merged.csv", 'w') as f:
# 	write = csv.writer(f)
# 	write.writerow(out_header)
# 	write.writerows(out)

# print("full file writing done")

# df = pd.read_csv("pNpS_merged.csv")
# df = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)]
# print(df[["gene_category", 'pn_ps']].groupby("gene_category").describe())
# print(df[["gene_category", 'pn_ps']].groupby("gene_category").mean(numeric_only=True))
# print(df[["gene_category", 'pn_ps']].groupby("gene_category").median(numeric_only=True))

