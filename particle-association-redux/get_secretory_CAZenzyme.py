import csv
import glob
from re import A

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

e_cutoff = float(1e-10)

all_CAZ_f = glob.glob("../../OM-RGC_v2_reference_pieces/dbcan_output_old/p_*/diamond.out")
secretory_CAZ = "./secretory_CAZenzyme.txt"
secretory_pep = "./secretory_peptidase.txt"

def read_CAZ(CAZ_f):
	print(CAZ_f)
	info = list(csv.reader(open(CAZ_f), delimiter='\t'))[1:]
	out = [] # [[Gene_ID, DB_code], [], ...]
	for l in info:
		if float(l[-2]) > e_cutoff:
			continue # filter out things more than e-val cutoff
		Gene_ID = l[0].split("_")[0]
		DB_code = l[1].split("|")[1] # only look at the first affiliation to determine whether it is secretory
		out.append([Gene_ID, DB_code])
	return out

def separate_secretory(clean_out, secretory_set):
	all, secretory = [], []
	for l in clean_out:
		Gene_ID, DB_code = l
		all.append(Gene_ID)
		if l[1] in secretory_set:
			secretory.append(DB_code)
	return all, secretory

def main():
	sect_CAZ = set(line.strip() for line in open(secretory_CAZ))
	# secr_pep = set(line.strip() for line in open(secretory_pep))

	root = '/researchdrive/zhongj2/BLAST_tara_OM-RGC_v2/'
	f1, f2 = 'all_CAZenzyme.txt', "secretory_CAZenzyme.txt"
	with open(root + f1, 'w') as a, open(root + f2, "w") as b:
		for f in all_CAZ_f:
			clean_CAZ = read_CAZ(f)
			CAZ_list, secretory_CAZ_list = separate_secretory(clean_CAZ, sect_CAZ)
			for Gene_ID in CAZ_list: a.write("%s\n" % Gene_ID)
			for Gene_ID in secretory_CAZ_list: b.write("%s\n" % Gene_ID)

	print('Done, file outputed to: ' + root)
	print('filenames: ' + f1 + " ; " + f2)

if __name__ == "__main__":
	main()

