# Written by Jimmy Zhong (zhongj2@carleton.edu), Carleton '23 under Professor Rika Anderson
# date: August 4th 2021
import glob
import os
import subprocess
import sys
import csv
import time
import re

def get_COG_num_from_dict(COG_dict, search_str):
	COG_numbers = set()
	for key in COG_dict:
		if search_str in key.lower():
			COG_numbers.add(COG_dict[key])
			print(key.lower())
	print(str(COG_numbers))
	return COG_numbers

def get_reverse_COG_dict(category_or_funcion = 3):
	''' get a COG function dict, 
		2 => COG category, 3 => function
		key = gene_categories (text), value = COG_number
	'''
	D = {}
	cogfile = open('/usr/local/data/Python_scripts/COG_numbers_categories_annotations.txt', 'r')
	for line in cogfile:
		columns = line.split('\t')
		cog_number = columns[0]
		cog_gene_name = columns[category_or_funcion]
		D[cog_gene_name] = cog_number
	cogfile.close()
	return D

def get_dict_bin_coverage(sam_console_output_file):
	file = open(sam_console_output_file, 'r')
	lines = file.readlines()

	bin_coverage_dict = {}

	for i in range(len(lines)):
		if "reads; of these:" in lines[i]:
			bin_name = (lines[i-1]).strip()

			total_reads = lines[i].replace("reads; of these:", "")
			total_reads = int(total_reads.strip())
			aligned_reads = lines[i+3].replace("%) aligned exactly 1 time", "")
			aligned_reads = int(aligned_reads[:-6]) # get rid of the number
			alignment_rate = aligned_reads/total_reads

			bin_coverage_dict[bin_name] = alignment_rate
	
	dict_items = bin_coverage_dict.items()
	first_ten = list(dict_items)[:10]

	print(first_ten)
	print(len(bin_coverage_dict))
	return bin_coverage_dict

src_directory = "/workspace/data/zhongj/Transposase_Project/deep_ocean_bins/acinas-et-al-2020_Malaspina-deep-HQ-MQ-317-MAG-set" # where the bins are
analysis_dir = "/workspace/data/zhongj/Transposase_Project/deep_ocean_bins/bins_analysis"
bin_function_dir = "/workspace/data/zhongj/Transposase_Project/deep_ocean_bins/acinas-et-al-2020_Malaspina-deep-HQ-MQ-317-MAG-set-annotation"
small_size_file = "all_bins_mapped_small_hole.txt"
big_size_file = "all_bins_mapped_big_hole.txt"

os.chdir(src_directory)
path = r'./'
list_done = []
list_to_do = glob.glob(os.path.join(path, '*.fasta'))

# COG_function_dict = get_reverse_COG_dict(3)
# COG_category_dict = get_reverse_COG_dict(2)
# TA_COG_num_set = get_COG_num_from_dict(COG_function_dict, "toxin")

output = subprocess.check_output("grep 'toxin' /usr/local/data/Python_scripts/COG_numbers_categories_annotations.txt | awk {'print $1'}", shell=True, stderr=subprocess.STDOUT)
TA_COG_num_set = set(output.decode('utf-8').splitlines())
print(TA_COG_num_set)

output = subprocess.check_output("grep 'Defense_mechanisms' /usr/local/data/Python_scripts/COG_numbers_categories_annotations.txt | awk {'print $1'}", shell=True, stderr=subprocess.STDOUT)
defense_COG_num_set = set(output.decode('utf-8').splitlines())
print(defense_COG_num_set)

os.chdir(analysis_dir)
small_size_bin_cov = get_dict_bin_coverage(small_size_file) 
big_size_bin_cov = get_dict_bin_coverage(big_size_file) 

def trace_bin_origin(bin_name):
	if bin_name in small_size_bin_cov:
		this_bin_cov_in_small = float(small_size_bin_cov[bin_name])
		this_bin_cov_in_big = float(big_size_bin_cov[bin_name])

		ratio = this_bin_cov_in_small/this_bin_cov_in_big
		if ratio > 10: 
			return "free_living", str(ratio)
		elif ratio < 0.1:
			return "particle-attached", str(ratio)
		else:
			return "undecided", str(ratio)
	return "error", str(-1)

def findUniqueORFs(fileName):
	matchedORFs = set()
	with open(fileName) as f:
		for line in f:
			sl = line.split('\t')
			matchedORFs.add(sl[1])
	# make into a list and sort
	# sorted_list = list(matchedORFs).sort()
	
	# I just need the length of the list
	return len(matchedORFs)

def get_biofilm_count(bin_name):
	os.chdir(analysis_dir + "/biofilm_blast")
	filename = f"Biofilm_vs_{bin_name}.blastp"
	return findUniqueORFs(filename)

def get_transposase_count(bin_name):
	os.chdir(analysis_dir + "/transposase_blast")
	filename = f"Transposase_vs_{bin_name}.blastp"
	return findUniqueORFs(filename)

def count_genes_in_cog_set(bin_name, cog_set):
	count = 0
	Total_gene_count = -1 # first
	os.chdir(bin_function_dir)
	tsv_file = open(f"./{bin_name}.enhanced.tsv", 'r')
	for line in tsv_file:
		Total_gene_count += 1
		columns = line.split('\t')
		cog_gene_name = columns[5]
		if cog_gene_name in cog_set:
			count += 1
	return count, Total_gene_count

final_write = ['''bin_name,size_fraction,small/big ratio,biofilm_count,transposase_count,toxin_antitoxin_count,defense_count,total_gene_count,biofilm_prop,trans_prop,TA_prop,defense_prop''']

for filename in list_to_do:
	MAG_file = filename.strip('./')
	MAG_name = MAG_file.replace('.fasta', '')
	MAG_origin, numeric_ratio = trace_bin_origin(MAG_name)
	TA_count, Total_gene_count = count_genes_in_cog_set(MAG_name, TA_COG_num_set)
	defense_count, Total_gene_count = count_genes_in_cog_set(MAG_name, defense_COG_num_set)
	Biofilm_count = get_biofilm_count(MAG_name)
	Trans_count = get_transposase_count(MAG_name)

	Biofilm_prop = Biofilm_count/Total_gene_count
	Trans_prop = Trans_count/Total_gene_count
	TA_prop = TA_count/Total_gene_count
	defense_prop = defense_count/Total_gene_count

	this_line = f'''{MAG_name},{MAG_origin},{numeric_ratio},{str(Biofilm_count)},{str(Trans_count)},{str(TA_count)},{str(defense_count)},{str(Total_gene_count)},{str(Biofilm_prop)},{str(Trans_prop)}, {str(TA_prop)}, {str(defense_prop)}'''
	
	final_write.append(this_line)

os.chdir(analysis_dir)
outfile = open('origin_biofilm_trans_TA_each_bin.csv', 'w')
for line in final_write:
	outfile.write(line + '\n')
outfile.close()