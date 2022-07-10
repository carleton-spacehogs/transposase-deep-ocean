import csv
import glob
from Bio import SeqIO

def clean_name(MAG_name):
	if "mp-deep_mag-" in MAG_name:
		return MAG_name.replace("mp-deep_mag-","deep_MAG_")
	else:
		return MAG_name.replace("TOBG_","")

def get_MAG_contig_length_dict(MAG_fna_f):
	contigs = SeqIO.parse(open(MAG_fna_f), 'fasta')
	seq_length_dict = { Mseq.id:len(Mseq.seq) for Mseq in contigs }
	return seq_length_dict

def cal_MAG_composition(euk_file, domains):
	composition_dict = {d:(0) for d in domains}
	MAG_name = euk_file.split("/")[-2]
	MAG_fna = euk_file.replace("/tiara_euk.txt",".fasta") if "deep" in euk_file else euk_file.replace("/tiara_euk.txt",".fna")
	MAG_contig_dict = get_MAG_contig_length_dict(MAG_fna)
	MAG_total_bp = sum(MAG_contig_dict.values())
	composition_info = list(csv.reader(open(euk_file), delimiter='\t'))[1:]
	for row in composition_info:
		contig_name, contig_domain = row[0], row[1].strip()
		composition_dict[contig_domain] += int(MAG_contig_dict[contig_name])
	prop_composition = [composition_dict[d]/MAG_total_bp for d in domains]
	out_row = [clean_name(MAG_name), MAG_total_bp] + prop_composition
	for d in domains:
		out_row.append(composition_dict[d])
	return out_row

def main():
	domains=["archaea","bacteria","eukarya","prokarya","organelle","unknown"]
	all_composition_files = glob.glob("../../bins/*/*/tiara_euk.txt")
	rows = []
	for f in all_composition_files:
		new_row = cal_MAG_composition(f,domains)
		print(new_row)
		rows.append(new_row)
	
	top_row=f"bin,total_bp,{','.join([d+'_prop' for d in domains])},{','.join(domains)}\n"
	rows = sorted(rows,key=lambda x: x[4])

	out_f = "MAG_composition_summary.csv"
	with open(out_f, 'w') as a:
		a.write(top_row)
		csv_write = csv.writer(a, delimiter=',')
		csv_write.writerows(rows)
	print(f"outfile: -- {out_f} --!!!")

if __name__ == "__main__":
	main()
