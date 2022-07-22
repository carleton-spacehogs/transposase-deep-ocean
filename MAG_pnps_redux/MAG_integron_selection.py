import pandas as pd
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, MAG_base_pnps, get_all_pnps, ocean_depths_with_deep

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

def get_integron_dict(integron_rt, ocean):
	# int_f =f"{integron_rt}/{ocean}.integrons"
	int_f =f"{integron_rt}/{ocean}/all_bins_db/Results_Integron_Finder_{ocean}_all_bins/all_integron.txt"
	intgron_f = pd.read_csv(int_f,comment="#",sep="\t")
	print(intgron_f)
	integrons = intgron_f[["ID_replicon","pos_beg","pos_end","type"]].dropna()
	integrons_2D = integrons.values.tolist()
	integron_dict = {}
	for l in integrons_2D:
		contig, start, stop, gene_type = l
		if start >= stop: 
			exit(f"error while parsing: {int_f} --- start later than stop")
		if contig in integron_dict:
			integron_dict[str(contig)].append([int(start), int(stop), str(gene_type)])
		else:
			integron_dict[str(contig)] = [[int(start), int(stop), str(gene_type)]]
	return integron_dict

def merge_anvi_gene_call(MAG_rt, ocean, integron_dict):
	MAG_f = f"{MAG_rt}/{ocean}/all_bins_db/gene_callers_id-contig_full.txt"
	MAGs = pd.read_csv(MAG_f, sep="\t")[["gene_callers_id","contig","start","stop"]]
	MAG_2D = MAGs.values.tolist()
	int_cols = ["gene_callers_id", "contig", "start", "stop", "integron", "ocean"]
	MAG_integrons = []
	for l in MAG_2D:
		gene_callers_id,contig,start,stop = l
		if start >= stop: 
			exit(f"error while parsing: {MAG_f} --- start later than stop")
		else:
			is_integron = determine_if_integron(integron_dict, contig, start, stop)
			if not is_integron == "non_integron":
				MAG_integrons.append([gene_callers_id, contig, start, stop, is_integron, ocean])
	MAG_integrons_df = pd.DataFrame(MAG_integrons, columns=int_cols)
	return MAG_integrons_df

def merge_pnps_all(df, o_depth, MAG_rt):
	all_pnps = pd.DataFrame()
	for ocean, depths in o_depth.items():
		ocean_pnps = get_all_pnps(depths, MAG_rt, ocean)
		ocean_pnps["ocean"] = ocean
		all_pnps = all_pnps.append(ocean_pnps)
	merged_integrons_str = df.astype(str)
	MAG_int_pnps = merged_integrons_str.merge(all_pnps, on=["gene_callers_id","ocean"])
	return MAG_int_pnps

def merge_COG_cate(df, o_depth, MAG_rt):
	all_COG_cate = pd.DataFrame()
	for ocean, depths in o_depth.items():
		COG_f=f"{MAG_rt}/{ocean}/all_bins_db/anvi_genes_COG_categories.tsv"
		COG_cate = pd.read_csv(COG_f, sep = "\t").astype(str)
		COG_cate["ocean"] = ocean
		all_COG_cate = all_COG_cate.append(COG_cate)
	df = df.merge(all_COG_cate, on=["gene_callers_id","ocean"])
	return df

# main():
o_depth = ocean_depths()
MAG_rt="/researchdrive/zhongj2/MAG_pnps_redux"
integron_rt="/researchdrive/zhongj2/integron_finder_tara_contig"

merged_integrons = pd.DataFrame()
for ocean in o_depth.keys():
	print(ocean)
	# integron_dict = get_integron_dict(integron_rt, ocean)
	integron_dict = get_integron_dict(MAG_rt, ocean)
	ocean_integrons = merge_anvi_gene_call(MAG_rt, ocean, integron_dict)
	merged_integrons = merged_integrons.append(ocean_integrons)
	print(ocean_integrons)

merged_integrons.to_csv(path_or_buf=f'anvi_integrons_all_oceans.csv', sep=',', index=False)
MAG_int_pnps = merge_pnps_all(merged_integrons, o_depth, MAG_rt)
MAG_int_pnps = merge_COG_cate(MAG_int_pnps, o_depth, MAG_rt)

bin_info = pd.read_csv("signalCAZyme-abun-vs-genes-pnps.csv").astype(str)
int_pnps_col = list(MAG_int_pnps.columns)
int_pnps_col.insert(1, "bin")

MAG_int_2D = MAG_int_pnps.values.tolist()
splited_init = [[l[0]] + l[1].split("_", 1) + l[2:] for l in MAG_int_2D]

MAG_int_pnps2 = pd.DataFrame(splited_init)
MAG_int_pnps2.columns = int_pnps_col

MAG_int_pnps2 = MAG_int_pnps2.merge(bin_info, on=["bin","sample_id"])

out = []
for l in MAG_int_pnps2.values.tolist():
	COG_category_pos = 11
	COG_cate = l[COG_category_pos]
	if "!!!" in COG_cate:
		all_COG_cates = COG_cate.split("!!!")
		for COG in all_COG_cates:
			out.append( l[:COG_category_pos] + [COG] + l[COG_category_pos+1:] )
	else:
		out.append(l)

MAG_int_pnps3 = pd.DataFrame(out)
MAG_int_pnps3.columns = list(MAG_int_pnps2.columns)

MAG_int_pnps3.to_csv(path_or_buf=f'MAG_integron_pnps.csv', sep=',', index=False)

# if __name__ == "__main__":
# 	main()


