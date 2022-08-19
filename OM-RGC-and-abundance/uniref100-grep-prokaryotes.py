from Bio import SeqIO
import csv

base="/researchdrive/zhongj2"
only_prok=f"{base}/uniref100_only_prok.fasta"

all_prok=list(csv.reader(open(f"{base}/NCBI_taxon/TaxID_category.txt"), delimiter=' '))
prok_dict = {l[0]:l[1] for l in all_prok} # taxon_id:domain

with open(only_prok, 'w') as out_file:
	fasta_seq = SeqIO.parse(open(f"{base}/uniref100.fasta"), 'fasta')
	for f in fasta_seq:
		des = f.description
		tax_id = des.split("TaxID=")[1].split(" ")[0]
		if tax_id in prok_dict:
			domain = prok_dict[tax_id]
			head = f">{domain}_{f.id}"
			out_file.write( f"{head}\n{str(f.seq)}\n" )

print(f"finish writing: {only_prok}")
