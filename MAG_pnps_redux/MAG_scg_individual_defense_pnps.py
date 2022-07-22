import sys
import pandas as pd
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, merge_COG_category, get_all_pnps, get_binning_info, merge_blastp_f, ocean_depths_with_deep

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

o_depths = ocean_depths()
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

defense_integron = individual_defense_pnps.merge(integrons, on=["ocean", "gene_callers_id"], how = "left")
d_concise = defense_integron[["gene_callers_id", "pnps", "ocean", "depth", "bin", "defense_count", 
"sample_id", "Trans_abun","scg_pnps_median","scg_count",'integron']]
d_concise.to_csv(path_or_buf=f'individual_defense_pnps.csv', sep=',', index=False)

# for transposase
individual_trans_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_trans = per_gene_ORFs_per_sample_pnps(ocean, root, depths, avail_MAG, "transposase")
	ocean_trans["ocean"] = ocean
	individual_trans_pnps = individual_trans_pnps.append(ocean_trans)

trans_concise = individual_trans_pnps[["gene_callers_id", "pnps", "ocean", "depth", "bin", "defense_count", 
"sample_id", "Trans_abun","scg_pnps_median","scg_count"]]
trans_concise.to_csv(path_or_buf=f'individual_transposase_pnps.csv', sep=',', index=False)

# # for toxin!!!
individual_toxin_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_toxin = per_gene_ORFs_per_sample_pnps(ocean, root, depths, avail_MAG, "toxin")
	ocean_toxin["ocean"] = ocean
	individual_toxin_pnps = individual_toxin_pnps.append(ocean_toxin)

t_concise = individual_toxin_pnps[["gene_callers_id", "pnps", "ocean", "depth", "bin", "defense_count", 
"sample_id", "signal_CAZyme_abun","scg_pnps_median","scg_count"]]
t_concise.to_csv(path_or_buf=f'individual_toxin_pnps.csv', sep=',', index=False)

