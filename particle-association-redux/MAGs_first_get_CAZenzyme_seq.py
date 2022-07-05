import csv
import glob
import sys
from first_get_CAZenzyme_and_peptidase import get_gene_amino_acid_base

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

def read_overview(CAZ_f):
	info = list(csv.reader(open(CAZ_f), delimiter='\t'))[1:]
	out = set()
	for l in info:
		if int(l[-1]) >= 2: # has have 2/3 tools agree
			out.add(l[0])
	return out

# def read_peptidase(pep_f): will develop later, probably import from first_get_CAZenzyme...py

def extract_specific_gene_sequences(overviews, aminoacids_fs):
	if len(overviews) != len(aminoacids_fs):
		uniInput_dir_set = set([a.replace("uniInput","") for a in aminoacids_fs])
		overview_dir_set = set([o.replace("overview.txt","") for o in overviews])
		print(uniInput_dir_set.difference(overview_dir_set))
		print(overview_dir_set.difference(uniInput_dir_set))
		exit("the length of the overviews and amino acids uniInput do not match")
	for i in range(len(overviews)):
		gene_id_set = read_overview(overviews[i])
		subset_aa = get_gene_amino_acid_base(aminoacids_fs[i], gene_id_set, "only_CAZenzyme.faa")
		print(f"write new file: {subset_aa}")

def main():
	ocean = sys.argv[1] # e.g. python3 MAGs_get_CAZenzyme_seq.py ARS
	path = f"../../bins/{ocean}/TOBG_{ocean}-*/"
	if ocean == "deep": path = "../../bins/deep/mp-deep_mag-*"

	all_CAZ_overview = glob.glob(path + "overview.txt")
	all_aminoacid_f = glob.glob(path + "uniInput")

	extract_specific_gene_sequences(all_CAZ_overview, all_aminoacid_f)
	print("done with CAZenzyme!!!")

if __name__ == "__main__":
	main()
