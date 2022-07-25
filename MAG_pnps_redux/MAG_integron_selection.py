import pandas as pd
import csv
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, MAG_base_pnps, get_all_pnps, get_binning_info

deep_oceans = ["IN","NAT","SAT","NP","SP"]

def within(start1, start2, end1, end2):
	if abs(int(start1)-int(start2)) < 100 and abs(int(end1)-int(end2)) < 100:
		return True
	else:
		return False

def determine_if_integron(integron_dict, MAG_contig, MAG_start, MAG_end):
	# cleaned_contig = MAG_contig.split('_', 1)[1]
	if MAG_contig in integron_dict:
		integron_list = integron_dict[MAG_contig]
		for integron in integron_list:
			beg, end, g_type = integron
			if within(MAG_start, beg, MAG_end, end):
				return g_type
	return "non_integron"

def integron_dict_helper(integron_rt, ocean, integron_dict):
	int_f =f"{integron_rt}/{ocean}/all_bins_db/Results_Integron_Finder_{ocean}_all_bins/all_integron.txt"
	intgron_f = pd.read_csv(int_f,comment="#",sep="\t")
	integrons = intgron_f[["ID_replicon","pos_beg","pos_end","type"]].dropna()
	integrons_2D = integrons.values.tolist()
	for l in integrons_2D:
		contig, start, stop, gene_type = l
		if start >= stop: 
			exit(f"error while parsing: {int_f} --- start later than stop")
		if contig in integron_dict:
			integron_dict[str(contig)].append([int(start), int(stop), str(gene_type)])
		else:
			integron_dict[str(contig)] = [[int(start), int(stop), str(gene_type)]]
	return integron_dict

def get_integron_dict(integron_rt, ocean):
	# int_f =f"{integron_rt}/{ocean}.integrons"
	ocean_integrons = {}
	ocean_integrons = integron_dict_helper(integron_rt, ocean, ocean_integrons)
	if ocean in deep_oceans:
		ocean_integrons = integron_dict_helper(integron_rt, "deep", ocean_integrons)
	return ocean_integrons

def full_binning_info_helper(fname, source):
	tmp = list(csv.reader(open(fname), delimiter="\t"))[1:]
	bin_info = [l[:4] + l[1].split('_', 1) for l in tmp]
	bin_info = pd.DataFrame(bin_info)
	bin_info.columns = ["gene_callers_id","contig","start","stop","bin","bin_contig"]
	bin_info['source'] = source
	return bin_info

def get_full_binning_info(root, ocean):
	tail="all_bins_db/gene_callers_id-contig_full.txt"
	bin_info = full_binning_info_helper(f"{root}/{ocean}/{tail}", 'tara')
	if ocean in deep_oceans:
		deep_bin = full_binning_info_helper(f"{root}/{ocean}/{tail}", 'malaspina')
		bin_info = bin_info.append(deep_bin)
	return bin_info

def merge_anvi_gene_call(MAG_rt, ocean, integron_dict):
	MAG = get_full_binning_info(MAG_rt, ocean)
	MAG_2D = MAG.values.tolist()
	int_cols = list(MAG.columns) + ["integron", "ocean"]
	MAG_integrons = []
	for l in MAG_2D:
		contig, start, stop = l[1:4]
		start, stop = int(start), int(stop)
		if start > stop:
			print(f"{contig}, {start}, {stop}")
			print(f"error while parsing: {ocean} --- start later than stop")
		else:
			is_integron = determine_if_integron(integron_dict, contig, start, stop)
			if not is_integron == "non_integron":
				MAG_integrons.append(l + [is_integron, ocean])
	MAG_integrons_df = pd.DataFrame(MAG_integrons, columns=int_cols)
	return MAG_integrons_df

def merge_pnps_all(df, o_depth, MAG_rt):
	all_pnps = pd.DataFrame()
	for ocean, depths in o_depth.items():
		ocean_pnps = get_all_pnps(depths, MAG_rt, ocean)
		ocean_pnps["ocean"] = ocean
		all_pnps = all_pnps.append(ocean_pnps)
	merged_integrons_str = df.astype(str)
	MAG_int_pnps = merged_integrons_str.merge(all_pnps, on=["gene_callers_id","ocean","source"])
	return MAG_int_pnps

def COG_cate_helper(MAG_rt, ocean):
	tail="all_bins_db/anvi_genes_COG_categories.tsv"
	COG_cate_df = pd.read_csv(f"{MAG_rt}/{ocean}/{tail}", sep = "\t")
	COG_cate_df['source'] = 'tara'
	if ocean in deep_oceans:
		deep_COG_cate = pd.read_csv(f"{MAG_rt}/deep/{tail}", sep = "\t")
		deep_COG_cate["source"] = "malaspina"
		COG_cate_df = COG_cate_df.append(deep_COG_cate)
	out = []
	for l in COG_cate_df.values.tolist():
		COG_category_pos = 2
		COG_cate = l[COG_category_pos]
		if "!!!" in COG_cate:
			all_COG_cates = COG_cate.split("!!!")
			for COG in all_COG_cates:
				out.append( l[:COG_category_pos] + [COG] + l[COG_category_pos+1:] )
		else:
			out.append(l)
	COG_out = pd.DataFrame(out)
	print(COG_out)
	print(COG_cate_df)
	COG_out.columns = list(COG_cate_df.columns)
	COG_out["ocean"] = ocean
	return COG_out.astype(str)

def merge_COG_cate(df, o_depth, MAG_rt):
	all_COG_cate = pd.DataFrame()
	for ocean, depths in o_depth.items():
		COG_cate = COG_cate_helper(MAG_rt, ocean)
		all_COG_cate = all_COG_cate.append(COG_cate)
	df = df.merge(all_COG_cate, on=["gene_callers_id","ocean","source"], how = "left")
	return df

# main():
o_depth = ocean_depths()
MAG_rt="/researchdrive/zhongj2/MAG_pnps_redux"
# integron_rt="/researchdrive/zhongj2/integron_finder_tara_contig"

merged_integrons = pd.DataFrame()
for ocean in o_depth.keys():
	print(ocean)
	integron_dict = get_integron_dict(MAG_rt, ocean)
	ocean_integrons = merge_anvi_gene_call(MAG_rt, ocean, integron_dict)
	merged_integrons = merged_integrons.append(ocean_integrons)

merged_integrons.to_csv(path_or_buf=f'anvi_integrons_all_oceans.csv', sep=',', index=False)
tmp = merge_pnps_all(merged_integrons, o_depth, MAG_rt)
MAG_int_pnps = merge_COG_cate(tmp, o_depth, MAG_rt)
MAG_int_pnps.to_csv(path_or_buf=f'MAG_integron_pnps_all.csv', sep=',', index=False)

bin_info = pd.read_csv("signalCAZyme-abun-vs-genes-pnps.csv").astype(str)
MAG_int_pnps2 = MAG_int_pnps.merge(bin_info, on=["bin","sample_id"])

MAG_int_pnps2.to_csv(path_or_buf=f'MAG_integron_pnps_with_high_cov_bins.csv', sep=',', index=False)

# if __name__ == "__main__":
# 	main()


