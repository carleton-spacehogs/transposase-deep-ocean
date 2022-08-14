#!/bin/python3
from first_get_CAZenzyme_and_peptidase import get_gene_amino_acid_base, read_peptidase
from MAGs_first_get_CAZyme_peptidase_seq import read_overview_2_agree

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux
def extract_CAZyme_sequences(ocean, root):
	overview = f"{root}/{ocean}_CAZyme/overview.txt"
	protein_f = f"{root}/{ocean}-anvi-gene-calls.faa"
	out_f = f"{root}/{ocean}-only_CAZyme.faa"
	gene_id_set = read_overview_2_agree(overview)
	subset_aa = get_gene_amino_acid_base(protein_f, gene_id_set, outfile = out_f, protein_fname = protein_f)
	print(subset_aa)

def extract_peptidase_sequences(ocean, root):
	protein_f = f"{root}/{ocean}/all_bins_db/anvi-gene-calls.faa"
	out_f = f"{root}/{ocean}/only_peptidase.faa"
	blastp = f"{root}/{ocean}/peptidase_diamond_unique.blastp"
	gene_id_set = read_peptidase(blastp, split=False)
	subset_aa = get_gene_amino_acid_base(protein_f, gene_id_set, outfile = out_f, protein_fname = protein_f)
	print(subset_aa)

def main():
	oceans = ["NAT","SAT","IN","NP","SP","ARS","CPC","EAC","MED","RS","deep"]
	root1="/workspace/data/zhongj/MAG_CAZymes"
	root2="/researchdrive/zhongj2/MAG_pnps_redux"
	for ocean in oceans:
		# extract_CAZyme_sequences(ocean, root1)
		extract_peptidase_sequences(ocean, root2)

if __name__ == "__main__":
	main()
