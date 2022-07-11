import csv
import glob
import sys
from first_get_CAZenzyme_and_peptidase import get_gene_amino_acid_base
from first_get_CAZenzyme_and_peptidase import read_peptidase

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

def check_sample_length(list1, list2, text1, text2):
	if len(list1) != len(list2):
		set_1 = set([a.replace(text1,"") for a in list1])
		set_2 = set([o.replace(text2,"") for o in list2])
		print(set_1.difference(set_2))
		print(set_2.difference(set_1))
		exit(f"the length of the {text1} and {text2} do not match")

def read_overview_only_diamond(CAZ_f):
	info = list(csv.reader(open(CAZ_f), delimiter='\t'))[1:]
	out = set()
	for l in info:
		if l[-2].strip() != "-": # diamond recognizes it -> only_diamond_CAZenzyme.faa
			out.add(l[0])
	return out

def read_overview_2_agree(CAZ_f):
	info = list(csv.reader(open(CAZ_f), delimiter='\t'))[1:]
	out = set()
	for l in info:
		if int(l[-1]) >= 2: # has have 2/3 tools agree -> only_CAZenzyme.faa
			out.add(l[0])
	return out

def extract_CAZenzyme_sequences(overviews, aminoacids_fs):
	check_sample_length(overviews, aminoacids_fs, "overview.txt", "uniInput")
	for i in range(len(overviews)):
		gene_id_set = read_overview_only_diamond(overviews[i])
		subset_aa = get_gene_amino_acid_base(aminoacids_fs[i], gene_id_set, "only_diamond_CAZenzyme.faa")
		print(f"write new file: {subset_aa}")

def extract_peptidase_sequences(blastps, aminoacids_fs):
	check_sample_length(blastps, aminoacids_fs, "peptidase_diamond_unique.blastp", "uniInput")
	for i in range(len(blastps)):
		gene_id_set = read_peptidase(blastps[i], split=False)
		subset_aa = get_gene_amino_acid_base(aminoacids_fs[i], gene_id_set, "only_diamond_peptidase.faa")
		print(f"write new file: {subset_aa}")

def main():
	ocean = sys.argv[1] # e.g. python3 MAGs_get_CAZenzyme_seq.py ARS
	path = f"../../bins/{ocean}/TOBG_{ocean}-*/"
	if ocean == "deep": path = "../../bins/deep/mp-deep_mag-*/"

	all_CAZ_overview = glob.glob(path + "overview.txt")
	all_pep_blastres = glob.glob(path + "peptidase_diamond_unique.blastp")
	all_aminoacid_f = glob.glob(path + "uniInput")
	print(len(all_aminoacid_f))

	extract_CAZenzyme_sequences(all_CAZ_overview, all_aminoacid_f)
	print("done with CAZenzyme!!!")
	extract_peptidase_sequences(all_pep_blastres, all_aminoacid_f)
	print("done with peptidase!!!")

if __name__ == "__main__":
	main()
