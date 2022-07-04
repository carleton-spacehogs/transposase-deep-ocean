import csv
import glob
from Bio import SeqIO

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

e_cutoff = float(1e-10)

all_CAZ_f = glob.glob("../../OM-RGC_v2_reference_pieces/dbcan_output_old/p_*/diamond.out")
all_pep_f = glob.glob("../../OM-RGC_v2_reference_pieces/blast_petidase/p_*_peptidases.tsv")
all_aminoacid_f = glob.glob("../../OM-RGC_v2_reference_pieces/dbcan_output_old/p_*/uniInput")

def read_CAZ(CAZ_f):
	info = list(csv.reader(open(CAZ_f), delimiter='\t'))[1:]
	out = set()
	for l in info:
		if float(l[-2]) > e_cutoff:
			continue # filter out things more than e-val cutoff
		Gene_ID = l[0].split("_")[0]
		out.add(Gene_ID)
	return out

def read_peptidase(pep_f):
	print(pep_f)
	info = list(csv.reader(open(pep_f), delimiter='\t'))
	out = set()
	for l in info:
		if float(l[-2]) > e_cutoff:
			continue # filter out things more than e-val cutoff
		Gene_ID = l[1].split("_")[0]
		out.add(Gene_ID)
	return out

def get_gene_amino_acid(in_files, gene_list_file, outfile_name):
	gene_id_set = set(line.strip() for line in open(gene_list_file))
	print(gene_list_file)
	for in_f in in_files:
		out_filename = in_f.replace("uniInput", outfile_name)
		print(out_filename)
		with open(out_filename, 'a') as out_file:
			fasta_seq = SeqIO.parse(open(in_f), 'fasta')
			for f in fasta_seq:
				name, seq = f.id.split("_")[0], str(f.seq)
				if name in gene_id_set:
					out_file.write(f">{name}\n")
					out_file.write(f"{seq}\n")
	print(f"write file {out_filename} done!")

def main():
	# root = '/researchdrive/zhongj2/BLAST_tara_OM-RGC_v2/'
	f1, f2 = './all_CAZenzyme.txt', './all_peptidase.txt'
	with open(f1, 'w') as a:
		for f in all_CAZ_f:
			clean_CAZ = read_CAZ(f)
			for Gene_ID in clean_CAZ:
				a.write("%s\n" % Gene_ID)
	
	print("done with CAZenzyme!!!")

	with open(f2, 'w') as b:
		for f in all_pep_f:
			clean_pep = read_peptidase(f)
			for Gene_ID in clean_pep:
				b.write("%s\n" % Gene_ID)

	print(f"Done, file outputed to: {f1} ; {f2}")

	# get_gene_amino_acid(all_aminoacid_f, f1, "only_CAZenzyme.faa")
	get_gene_amino_acid(all_aminoacid_f, f2, "only_peptidase.faa")

if __name__ == "__main__":
	main()
