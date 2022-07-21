import sys
import pandas as pd
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, merge_COG_category, get_all_pnps, get_binning_info

def per_defense_ORFs_per_sample_pnps(ocean, root, depths):
	bin_info = get_binning_info(root, ocean)
	all_pnps = get_all_pnps(depths, root, ocean)
	all_pnps = merge_COG_category(all_pnps, root, ocean)
	defense_pnps = all_pnps[all_pnps['category_accession'].str.contains("V", na=False)]
	defense_pnps = defense_pnps.merge(bin_info, on = "gene_callers_id")
	defense_pnps = defense_pnps.merge(avail_MAG, how = "left", on = ["bin", "sample_id"]).reindex()
	return defense_pnps

o_depths = ocean_depths()
root="/researchdrive/zhongj2/MAG_pnps_redux"

all_MAG_pnps = pd.read_csv("signalCAZyme-abun-vs-genes-pnps.csv")
avail_MAG = all_MAG_pnps[all_MAG_pnps.scg_count >= 11].astype(str).drop("depth", axis = 1)
integrons = pd.read_csv("anvi_integrons_all_oceans.csv").astype(str)

individual_defense_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_defense = per_defense_ORFs_per_sample_pnps(ocean, root, depths)
	ocean_defense["ocean"] = ocean
	individual_defense_pnps = individual_defense_pnps.append(ocean_defense)

defense_integron = individual_defense_pnps.merge(integrons, on=["ocean", "gene_callers_id"], how = "left")

concise = defense_integron[["gene_callers_id", "pnps", "ocean", "depth", "bin", "defense_count",
"sample_id", "signal_CAZyme_abun","scg_pnps_median","ocean",'integron']]

individual_defense_pnps.to_csv(path_or_buf=f'individual_defense_pnps.csv', sep=',', index=False)

