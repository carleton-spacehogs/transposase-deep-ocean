import sys
import pandas as pd
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, merge_COG_category, get_all_pnps, get_binning_info, merge_blastp_f, read_pnps
from utils import MAG_db_fun

def get_all_pnps_v2(depths, root, ocean):
	pnps_all_sam = pd.DataFrame()
	if ocean == "deep":
		regions = ["IN", "NP", "NAT", "SAT", "SP"]
		for reg in regions:
			
	else:
		for depth in depths:
			tara_MAGs_pnps = read_pnps(root, ocean, depth)
			tara_MAGs_pnps['depth'] = str(depth)
			tara_MAGs_pnps['source'] = 'tara'
			pnps_all_sam = pnps_all_sam.append(tara_MAGs_pnps)
	pnps_all_sam.sort_values(by=['gene_callers_id'])
	return pnps_all_sam

def per_defense_ORFs_per_sample_pnps(ocean, root, depths, avail_MAG):
	bin_info = get_binning_info(root, ocean)
	all_pnps = get_all_pnps(depths, root, ocean)
	all_pnps = merge_COG_category(all_pnps, root, ocean)
	defense_pnps = all_pnps[all_pnps['category_accession'].str.contains("V", na=False)]
	defense_pnps = defense_pnps.merge(bin_info, on = ["gene_callers_id", 'source'])
	defense_pnps = defense_pnps.merge(avail_MAG, how = "inner", on = ["bin", "sample_id"]).reindex()
	return defense_pnps

def per_gene_ORFs_per_sample_pnps(ocean, root, depths, avail_MAG, gene):
	bin_info = get_binning_info(root, ocean)
	all_pnps = get_all_pnps(depths, root, ocean)
	toxin_pnps = merge_blastp_f(all_pnps, root, ocean, gene, merge_method = "inner")
	toxin_pnps = toxin_pnps.merge(bin_info, on = ["gene_callers_id", 'source'])
	toxin_pnps = toxin_pnps.merge(avail_MAG, how = "inner", on = ["bin", "sample_id"]).reindex()
	return toxin_pnps

o_depths = MAG_db_fun()
root="/researchdrive/zhongj2/MAG_pnps_redux"

# all_MAG_pnps = pd.read_csv("signalCAZyme-abun-vs-genes-pnps.csv")
all_MAG_pnps = pd.read_csv("transposase-abun-vs-scg-pnps.csv")
avail_MAG = all_MAG_pnps[all_MAG_pnps.scg_count >= 11].astype(str).drop("depth", axis = 1)
integrons = pd.read_csv("anvi_integrons_all_oceans.csv").astype(str)

# for the defense mechansims
individual_defense_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_defense = per_defense_ORFs_per_sample_pnps(ocean, root, depths, avail_MAG)
	ocean_defense["ocean"] = ocean
	individual_defense_pnps = individual_defense_pnps.append(ocean_defense)

id_col = ["ocean", "gene_callers_id"]
integron_tmp = integrons[ id_col + ["integron"] ]

out_col = ["gene_callers_id", "pnps", "ocean", "depth", "bin", "defense_count","sample_id",
"Trans_abun","scg_pnps_median","scg_count","MAG_pnps_median","total_count"]

defense_integron = individual_defense_pnps.merge(integron_tmp, on=id_col, how = "left")
d_concise = defense_integron[ out_col + ['integron'] ]
d_concise.to_csv(path_or_buf=f'individual_defense_pnps.csv', sep=',', index=False)

# for transposase
individual_trans_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_trans = per_gene_ORFs_per_sample_pnps(ocean, root, depths, avail_MAG, "transposase")
	ocean_trans["ocean"] = ocean
	individual_trans_pnps = individual_trans_pnps.append(ocean_trans)

trans_concise = individual_trans_pnps[out_col]
trans_concise.to_csv(path_or_buf=f'individual_transposase_pnps.csv', sep=',', index=False)

# # for toxin!!!
# individual_toxin_pnps = pd.DataFrame()
# for ocean, depths in o_depths.items():
# 	print(ocean)
# 	ocean_toxin = per_gene_ORFs_per_sample_pnps(ocean, root, depths, avail_MAG, "toxin")
# 	ocean_toxin["ocean"] = ocean
# 	individual_toxin_pnps = individual_toxin_pnps.append(ocean_toxin)

# t_concise = individual_toxin_pnps[["gene_callers_id", "pnps", "ocean", "depth", "bin", "defense_count", 
# "sample_id", "Trans_abun","scg_pnps_median","scg_count"]]
# t_concise.to_csv(path_or_buf=f'individual_toxin_pnps.csv', sep=',', index=False)

