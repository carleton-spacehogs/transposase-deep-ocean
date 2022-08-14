import sys
import pandas as pd
from MAG_singleCopyGene_defense_pnps_ratio import merge_COG_category, get_binning_info
from utils import get_bin_helper, read_pnps, MAG_db_fun, has_deep, data_root, merge_blastp_f

def get_all_pnps_v2(depths, root, ocean):
	pnps_all_sam = pd.DataFrame()
	if ocean == "deep":
		for reg in has_deep:
			mala_MAGs_pnps = read_pnps(root, ocean, reg)
			mala_MAGs_pnps['depth'] = 'deep'
			mala_MAGs_pnps['source'] = 'malaspina'
			pnps_all_sam = pnps_all_sam.append(mala_MAGs_pnps)
	else:
		for depth in depths:
			tara_MAGs_pnps = read_pnps(root, ocean, depth)
			tara_MAGs_pnps['depth'] = str(depth)
			tara_MAGs_pnps['source'] = 'tara'
			pnps_all_sam = pnps_all_sam.append(tara_MAGs_pnps)
	pnps_all_sam.sort_values(by=['gene_callers_id'])
	return pnps_all_sam

def per_defense_ORFs_per_sample_pnps(ocean, root, depths, MAG_db):
	bin_info = get_bin_helper(f"{root}/{ocean}/all_bins_db/gene_callers_id-contig.txt")
	all_pnps = get_all_pnps_v2(depths, root, ocean)
	all_pnps = merge_COG_category(all_pnps, root, ocean)
	defense_pnps = all_pnps[all_pnps['category_accession'].str.contains("V", na=False)]
	defense_pnps = defense_pnps.merge(bin_info, on="gene_callers_id")
	defense_pnps = defense_pnps.merge(MAG_db, how = "inner", on = ["bin", "sample_id", "depth"]).reindex()
	return defense_pnps

def per_gene_ORFs_per_sample_pnps(ocean, root, depths, MAG_db, gene):
	bin_info = get_bin_helper(f"{root}/{ocean}/all_bins_db/gene_callers_id-contig.txt")
	all_pnps = get_all_pnps_v2(depths, root, ocean)
	gene_pnps = merge_blastp_f(all_pnps, root, ocean, gene, merge_method = "inner")
	gene_pnps = gene_pnps.merge(bin_info, on = "gene_callers_id")
	gene_pnps = gene_pnps.merge(MAG_db, how = "inner", on = ["bin", "sample_id", "depth"]).reindex()
	return gene_pnps

def add_integrons(id_col, df):
	integrons = pd.read_csv("../integron_finder_v2/MAG-integron-all.csv").astype(str)
	integron_tmp = integrons[ id_col + ["integron"] ]
	df_integron = df.merge(integron_tmp, on=id_col, how = "left")
	return df_integron

o_depths = MAG_db_fun()
root="/researchdrive/zhongj2/MAG_pnps_redux"

all_MAG_pnps = pd.read_csv("MAG_info_db/MAG_all.csv")

# for the defense mechansims
individual_defense_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_defense = per_defense_ORFs_per_sample_pnps(ocean, data_root, depths, all_MAG_pnps)
	ocean_defense["ocean"] = ocean
	individual_defense_pnps = individual_defense_pnps.append(ocean_defense)

id_col = ["ocean", "gene_callers_id"]
out_col = ["gene_callers_id", "pnps", "ocean", "depth", "bin","sample_id",
"scg_pnps_median","scg_count","MAG_pnps_median","total_count", 'integron']

defense_integron = add_integrons(id_col, individual_defense_pnps)
defense_integron[out_col].to_csv(path_or_buf=f'individual_defense_pnps.csv', sep=',', index=False)

# for transposase
individual_trans_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_trans = per_gene_ORFs_per_sample_pnps(ocean, root, depths, all_MAG_pnps, "transposase")
	ocean_trans["ocean"] = ocean
	individual_trans_pnps = individual_trans_pnps.append(ocean_trans)

trans_integron = add_integrons(id_col, individual_trans_pnps)
trans_integron[out_col].to_csv(path_or_buf=f'individual_transposase_pnps.csv', sep=',', index=False)

# for toxin!!!
individual_toxin_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_toxin = per_gene_ORFs_per_sample_pnps(ocean, root, depths, all_MAG_pnps, "toxin")
	ocean_toxin["ocean"] = ocean
	individual_toxin_pnps = individual_toxin_pnps.append(ocean_toxin)

toxin_integron = add_integrons(id_col, individual_toxin_pnps)
toxin_integron[out_col].to_csv(path_or_buf=f'individual_toxin_pnps.csv', sep=',', index=False)
