import sys
import pandas as pd
from MAG_singleCopyGene_defense_pnps_ratio import merge_COG_category
import utils

def add_integrons(df):
	integrons = pd.read_csv("../integron_finder_v2/MAG-integron-all.csv").astype(str)
	integrons['source'] = "malaspina"
	integrons.loc[integrons['ocean'] != 'deep', 'source'] = 'tara-' + integrons['ocean']
	integron_tmp = integrons[["source", "gene_callers_id", "integron"]]
	df_integron = df.merge(integron_tmp, on=["source", "gene_callers_id"], how = "left")
	return df_integron

def get_all_pnps_v2(depths, root, ocean):
	pnps_all_sam = pd.DataFrame()
	if ocean == "deep":
		for reg in utils.has_deep:
			mala_MAGs_pnps = utils.read_pnps(root, ocean, reg)
			mala_MAGs_pnps[['depth', 'ocean', 'source']] = ['deep', reg, 'malaspina']
			pnps_all_sam = pnps_all_sam.append(mala_MAGs_pnps)
	else:
		for depth in depths:
			tara_MAGs_pnps = utils.read_pnps(root, ocean, depth)
			tara_MAGs_pnps[['depth', 'ocean', 'source']] = [depth, ocean, 'tara-' + ocean]
			pnps_all_sam = pnps_all_sam.append(tara_MAGs_pnps)
	pnps_all_sam.sort_values(by=['gene_callers_id'])
	return pnps_all_sam

# V: defense mechanism; T: signal transduction
def per_category_ORFs_per_sample_pnps(ocean, root, depths, MAG_db, category):
	bin_info = utils.get_bin_helper(f"{root}/{ocean}/all_bins_db/gene_callers_id-contig.txt")
	all_pnps = get_all_pnps_v2(depths, root, ocean)
	all_pnps = merge_COG_category(all_pnps, root, ocean)
	defense_pnps = all_pnps[all_pnps['category_accession'].str.contains(category, na=False)]
	defense_pnps = defense_pnps.merge(bin_info, on="gene_callers_id")
	defense_pnps = defense_pnps.merge(MAG_db, how = "inner", on = ["bin", "sample_id", "depth"]).reindex()
	return defense_pnps

# gene: toxin, transposase
def per_gene_ORFs_per_sample_pnps(ocean, root, depths, MAG_db, gene):
	bin_info = utils.get_bin_helper(f"{root}/{ocean}/all_bins_db/gene_callers_id-contig.txt")
	all_pnps = get_all_pnps_v2(depths, root, ocean)
	gene_pnps = utils.merge_blastp_f(all_pnps, root, ocean, gene, merge_method = "inner")
	gene_pnps = gene_pnps.merge(bin_info, on = "gene_callers_id")
	gene_pnps = gene_pnps.merge(MAG_db, how = "inner", on = ["bin", "sample_id", "depth"]).reindex()
	return gene_pnps

def secreting_ORFs_per_sample_pnps(ocean, root, depths, MAG_db):
	bin_info = utils.get_bin_helper(f"{root}/{ocean}/all_bins_db/gene_callers_id-contig.txt")
	secret_CAZ_df = utils.signalp_to_df_helper(f"{utils.signal_root}/{ocean}-CAZenzyme_gramPos_summary.signalp5")
	secret_info = secret_CAZ_df.merge(bin_info, on = "gene_callers_id")
	all_pnps = get_all_pnps_v2(depths, root, ocean)
	secret_pnps = secret_info.merge(all_pnps, on="gene_callers_id")
	secret_pnps = secret_pnps.merge(MAG_db, how = "inner", on = ["bin", "sample_id", "depth"]).reindex()
	return secret_pnps

o_depths = utils.MAG_db_fun()

id_col = ["ocean", "gene_callers_id"]
out_col = ["gene_callers_id", "pnps", "ocean", "depth", "bin","sample_id",
"scg_pnps_median","scg_count","MAG_pnps_median","total_count", 'integron']

all_MAG_pnps = pd.read_csv("MAG_info_db/MAG_all.csv")

# for the defense mechansims
individual_defense_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_defense = per_category_ORFs_per_sample_pnps(ocean, utils.data_root, depths, all_MAG_pnps, "V")
	individual_defense_pnps = individual_defense_pnps.append(ocean_defense)

defense_integron = add_integrons(individual_defense_pnps)
defense_integron[out_col].to_csv(path_or_buf=f'individual_defense_pnps.csv', sep=',', index=False)

# for the signal tranductions
individual_signal_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_signal = per_category_ORFs_per_sample_pnps(ocean, utils.data_root, depths, all_MAG_pnps, "T")
	individual_signal_pnps = individual_signal_pnps.append(ocean_signal)

signal_integron = add_integrons(individual_signal_pnps)
signal_integron[out_col].to_csv(path_or_buf=f'individual_signal_pnps.csv', sep=',', index=False)

# for the sectory genes (biofilms)
individual_secret_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_secret = secreting_ORFs_per_sample_pnps(ocean, utils.data_root, depths, all_MAG_pnps)
	individual_secret_pnps = individual_secret_pnps.append(ocean_secret)

secret_integron = add_integrons(individual_secret_pnps)
secret_integron[out_col].to_csv(path_or_buf=f'individual_secretary_pnps.csv', sep=',', index=False)

# for transposase
individual_trans_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_trans = per_gene_ORFs_per_sample_pnps(ocean, utils.data_root, depths, all_MAG_pnps, "transposase")
	individual_trans_pnps = individual_trans_pnps.append(ocean_trans)

trans_integron = add_integrons(individual_trans_pnps)
trans_integron[out_col].to_csv(path_or_buf=f'individual_transposase_pnps.csv', sep=',', index=False)

# for toxin!!!
individual_toxin_pnps = pd.DataFrame()
for ocean, depths in o_depths.items():
	print(ocean)
	ocean_toxin = per_gene_ORFs_per_sample_pnps(ocean, utils.data_root, depths, all_MAG_pnps, "toxin")
	individual_toxin_pnps = individual_toxin_pnps.append(ocean_toxin)

toxin_integron = add_integrons(individual_toxin_pnps)
toxin_integron[out_col].to_csv(path_or_buf=f'individual_toxin_pnps.csv', sep=',', index=False)
