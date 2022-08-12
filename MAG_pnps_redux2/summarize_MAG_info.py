import csv
import pandas as pd
import numpy as np
from utils import MAG_db_fun, data_root

MAG_db_col = ["bin", "sample_id", "depth", "source", "scg_pnps_median","scg_count","MAG_pnps_median","total_count"]

def get_bin_helper(fname):
	tmp = list(csv.reader(open(fname), delimiter=" "))[1:]
	bin_info = [[l[0]] + l[1].split('_', 1) for l in tmp]
	bin_info = pd.DataFrame(bin_info)
	bin_info.columns = ["gene_callers_id", "bin", "contig"]
	return bin_info

def get_binning_info(root, ocean, delim = " "):
	bin_info = get_bin_helper(f"{root}/{ocean}/all_bins_db/gene_callers_id-contig.txt")
	bin_info['source'] = 'tara' if ocean != "deep" else "malaspina"
	return bin_info

# gene = "toxin", "transposase"
def merge_blastp_f(df, root, ocean, gene, merge_method = "left"):
	gene_blastp=pd.read_csv(f"{root}/{ocean}/{gene}_diamond_unique.blastp", sep = "\t", header=None).astype(str)[[1,0]]
	gene_blastp.rename(columns={1:"gene_callers_id", 0:f"{gene}_id"}, inplace=True)
	df = df.merge(gene_blastp, on="gene_callers_id", how = merge_method)
	return df

def merge_single_copy_gens(df, root, ocean):
	base=f"{root}/{ocean}/all_bins_db"
	scg = pd.read_csv(f"{base}/hmm-call-Bacteria_71.txt", sep = "\t").astype(str)
	scg.rename(columns={"gene_name":"scg"}, inplace=True)
	scg = scg[["gene_callers_id","scg"]]
	out = df.merge(scg, on="gene_callers_id", how = "left")
	return out

def read_pnps(root, ocean, depth): # ocean and depth is reversed for the Malaspina deep ocean
	cols = ["gene_callers_id", "pnps", "sample_id"]
	depth_pnps = pd.read_csv(f"{root}/{ocean}/PROFILE-{depth}/all_MAGs_pnps/pNpS.txt", sep = "\t").astype(str)
	depth_pnps = depth_pnps[depth_pnps['pNpS_gene_reference'].astype(float) < 5] # sanity check
	depth_pnps.rename(columns= {"corresponding_gene_call":cols[0],"pNpS_gene_reference":"pnps"}, inplace=True)
	return depth_pnps[cols]

def get_all_pnps(depths, root, ocean):
	pnps_all_sam = pd.DataFrame()
	if ocean == "deep":
		for region in ["IN", "SAT", "NAT", "SP", "NP"]:
			deep_MAGs_pnps = read_pnps(root, ocean, depth = region)
			pnps_all_sam = pnps_all_sam.append(deep_MAGs_pnps)
		pnps_all_sam[['ocean','sample_id']] = pnps_all_sam['sample_id'].str.split('_',expand=True)
		pnps_all_sam['source'] = 'malaspina'
		pnps_all_sam['depth'] = 'deep'
	else:
		for depth in depths:
			tara_MAGs_pnps = read_pnps(root, ocean, depth)
			tara_MAGs_pnps['depth'] = str(depth)
			pnps_all_sam = pnps_all_sam.append(tara_MAGs_pnps)
		pnps_all_sam['source'] = 'tara'
	pnps_all_sam.sort_values(by=['gene_callers_id'])
	return pnps_all_sam

# def get_pnps_summary(pnps_df, colname, gene_name):
# 	gene_pnps = pnps_df[~pnps_df[colname].isnull()]
# 	summary = gene_pnps[["bin","sample_id","pnps"]].groupby(['bin','sample_id']).agg(
# 		median=('pnps', np.nanmedian),count=('pnps', np.size)).reset_index()
# 	summary.rename(columns={"median":f"{gene_name}_pnps_median","count":f"{gene_name}_count"}, inplace=True)
# 	return summary

def MAG_base_pnps(root, ocean, depths):
	bin_base = get_binning_info(root, ocean)
	bin_base = merge_single_copy_gens(bin_base, root, ocean)
	bin_base = merge_blastp_f(bin_base, root, ocean, "transposase")
	all_pnps = get_all_pnps(depths, root, ocean)
	pnps = all_pnps.merge(bin_base, on = ["gene_callers_id","source"], how = "left")
	pnps.pnps = pnps.pnps.astype(float)
	whole_MAG_pnps = pnps.groupby(['bin','sample_id']).agg(
		MAG_pnps_median=('pnps', np.nanmedian),
		total_count=('pnps', np.size)
	).reset_index()
	scg_summary = pnps[~pnps.scg.isnull()].groupby(['bin','sample_id']).agg(
		scg_pnps_median=('pnps', np.nanmedian),
		scg_count=('pnps', np.size)
	).reset_index()
	out = whole_MAG_pnps.merge(scg_summary, on=["bin","sample_id"], how = "left")
	return out

def main():
	MAG_db = MAG_db_fun()
	all_MAG = pd.DataFrame()
	final_out = "MAG_info_db/MAG_all.csv"
	for ocean, depths in MAG_db.items():
		ocean_db = MAG_base_pnps(data_root, ocean, depths)
		ocean_db.to_csv(f"MAG_info_db/MAG_{ocean}.csv", index=False)
		all_MAG = all_MAG.append(ocean_db)

	all_MAG.to_csv(final_out, index=False)
	print(f"wrote -------- {final_out} -------- !!!")

if __name__ == "__main__":
	main()

