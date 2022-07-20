import pandas as pd
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, MAG_base_pnps, get_all_pnps

def within(start1, start2, end1, end2):
	if abs(int(start1)-int(start2)) < 100 and abs(int(end1)-int(end2)) < 100:
		return True
	else:
		return False

def determine_if_integron(integron_dict, MAG_contig, MAG_start, MAG_end):
	cleaned_contig = MAG_contig.split('_', 1)[1]
	if cleaned_contig in integron_dict:
		beg, end, g_type = integron_dict[cleaned_contig]
		if within(MAG_start, beg, MAG_end, end):
			return g_type
	return "non_integron"

def get_integron_dict(integron_rt, ocean):
	int_f =f"{integron_rt}/{ocean}.integrons"
	intgron_f = pd.read_csv(int_f,comment="#",sep="\t")
	print(intgron_f)
	integrons = intgron_f[["ID_replicon","pos_beg","pos_end","type"]].dropna()
	integrons_2D = integrons.values.tolist()
	integron_dict = {}
	for l in integrons_2D:
		contig, start, stop, gene_type = l
		if start >= stop: 
			exit(f"error while parsing: {int_f} --- start later than stop")
		else:
			integron_dict[str(contig)] = [int(start), int(stop), str(gene_type)]
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

o_depth = ocean_depths()
MAG_rt="/researchdrive/zhongj2/MAG_pnps_redux"
integron_rt="/researchdrive/zhongj2/integron_finder_tara_contig"

# ocean="IN"
merged_integrons = pd.DataFrame()
for ocean in o_depth.keys():
	print(ocean)
	integron_dict = get_integron_dict(integron_rt, ocean)
	ocean_integrons = merge_anvi_gene_call(MAG_rt, ocean, integron_dict)
	merged_integrons = merged_integrons.append(ocean_integrons)

merged_integrons.to_csv(path_or_buf=f'anvi_integrons_all_oceans.csv', sep=',', index=False)

all_pnps = pd.DataFrame()
for ocean, depths in o_depth.items():
	ocean_pnps = get_all_pnps(depths, MAG_rt, ocean)
	ocean_pnps["ocean"] = ocean
	all_pnps = all_pnps.append(ocean_pnps)

all_pnps_str = all_pnps.astype(str)

merged_integrons_str = merged_integrons.astype(str)

MAG_int_pnps = merged_integrons_str.merge(all_pnps, on=["gene_callers_id","ocean"])

# too liitle samples to work with...
# failed. done
MAG_int_pnps.to_csv(path_or_buf=f'MAG_integron_pnps.csv', sep=',', index=False)



