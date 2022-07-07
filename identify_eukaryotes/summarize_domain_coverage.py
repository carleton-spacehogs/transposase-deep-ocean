import csv
import glob

# def cal_tara_composition(composition_info, euk_file_name, domain):
# 	composition_dict = {d:float(0) for d in domain}
# 	for row in composition_info:
# 		if len(row) < 3: continue
# 		domain = row[1].strip()
# 		contig_name_list = row[0].split("_")
# 		length, avg_cov = float(contig_name_list[-3]), float(contig_name_list[-1])
# 		composition_dict[domain] += length*avg_cov # total coverage in base pairs
# 	return composition_dict

def get_deep_per_contig_cov(samtool_idxstats_f):
	file = list(csv.reader(open(samtool_idxstats_f), delimiter='\t'))
	cov_dict={}
	for line in file:
		if len(line) > 3:
			contig_name, contig_len, num_reads, zero = line
			cov_dict[contig_name.strip()] = int(num_reads)
	return cov_dict

def cal_deep_composition(composition_info, euk_file_name, domain):
	composition_dict = {d:float(0) for d in domain}
	SRR_num = euk_file_name.split("_")[0] if "SRR" in euk_file_name else euk_file_name.split("_")[1]
	contig_cov_file = f"get_contig_cov/{SRR_num}_contig_cov.txt"
	contig_cov_dict = get_deep_per_contig_cov(contig_cov_file)
	for row in composition_info:
		contig_name, domain = row[0], row[1].strip()
		composition_dict[domain] += float(contig_cov_dict[contig_name])
	return composition_dict

def summarize_composition(euk_file, tara_or_deep_func, domain):
	composition_info = list(csv.reader(open(euk_file), delimiter='\t'))[1:]
	composition_dict = tara_or_deep_func(composition_info, euk_file, domain)
	dict_sum = sum(composition_dict.values())
	depth = "BAT"
	size_fraction = "0.2-0.8"
	if "Tara" in euk_file:
		depth = euk_file.split("_")[-4]
		size_fraction = euk_file.split("_")[-3]
	if "size08" in euk_file:
		size_fraction = "0.8-0.5"
	out = [euk_file.replace("_euk.txt", ""), depth, size_fraction, dict_sum]
	for d in domain:
		out.append(composition_dict[d]/dict_sum)
	return out

def main():
	domains=["archaea","bacteria","eukarya","prokarya","organelle","unknown"]
	all_composition_files = glob.glob("*_euk.txt")
	rows = []
	for f in all_composition_files:
		print(f)
		funct_in = cal_deep_composition # cal_tara_composition() is deprecated -> use mapping based approach instead
		row = summarize_composition(f, funct_in, domains)
		rows.append(row)
	
	rows = sorted(rows,key=lambda x: x[6])

	sum_f_col = "sample,depth,size_fraction,sum_cov_bp,archaea,bacteria,eukarya,prokarya,organelle,unknown\n"
	out_f = "domain_composition_summary2.csv"
	with open(out_f, 'w') as a:
		a.write(sum_f_col)
		csv_write = csv.writer(a, delimiter=',')
		csv_write.writerows(rows)
	print(f"outfile: -- {out_f} --!!!")

if __name__ == "__main__":
	main()
