#!/bin/python3
import csv
import pandas as pd
import sys
import numpy as np

def get_binning_info(root, ocean):
	tmp = list(csv.reader(open(f"{root}/{ocean}/all_bins_db/gene_callers_id-contig.txt"), delimiter=' '))[1:]
	bin_info = [[l[0]] + l[1].split('_', 1) for l in tmp]
	bin_info = pd.DataFrame(bin_info)
	bin_info.columns = ["gene_callers_id", "bin", "contig"]
	return bin_info

def merge_toxin(df, root, ocean):
	toxin=pd.read_csv(f"{root}/{ocean}/toxin_diamond_unique.blastp", sep = "\t", header=None).astype(str)[[1,0]]
	toxin.rename(columns={1:"gene_callers_id", 0:"toxin_id"}, inplace=True)
	df = df.merge(toxin, on="gene_callers_id", how = "left")
	return df

def merge_COG_category(df, root, ocean):
	base=f"{root}/{ocean}/all_bins_db"
	COG_category = pd.read_csv(f"{base}/anvi_genes_COG_categories.tsv", sep = "\t").astype(str)
	# print(COG_category)
	# print(df)
	df = df.merge(COG_category, on="gene_callers_id", how = "left")
	return df

def merge_single_copy_gens(df, root, ocean):
	base=f"{root}/{ocean}/all_bins_db"
	scg = pd.read_csv(f"{base}/hmm-call-Bacteria_71.txt", sep = "\t").astype(str)
	scg.rename(columns={"gene_name":"scg"}, inplace=True)
	scg = scg[["gene_callers_id","scg"]]
	out = df.merge(scg, on="gene_callers_id", how = "left")
	return out

def get_all_pnps(depths, root, ocean):
	cols = ["gene_callers_id", "pnps", "sample_id","depth"]
	pnps_all_sam = pd.DataFrame()
	for depth in depths:
		depth_pnps = pd.read_csv(f"{root}/{ocean}/PROFILE-{depth}/all_MAGs_pnps/pNpS.txt", sep = "\t").astype(str)
		depth_pnps = depth_pnps[depth_pnps['pNpS_gene_reference'].astype(float) < 5] # sanity check
		depth_pnps['depth'] = str(depth)
		depth_pnps.rename(columns= {"corresponding_gene_call":cols[0],"pNpS_gene_reference":"pnps"}, inplace=True)
		pnps_all_sam = pnps_all_sam.append(depth_pnps[cols])
	pnps_all_sam.sort_values(by=['gene_callers_id'])
	return pnps_all_sam

def get_pnps_summary(pnps_df, colname, gene_name):
	gene_pnps = pnps_df[~pnps_df[colname].isnull()]
	summary = gene_pnps[["bin","sample_id","pnps"]].groupby(['bin','sample_id']).agg(
		median=('pnps', np.nanmedian),count=('pnps', np.size)).reset_index()
	summary.rename(columns={"median":f"{gene_name}_pnps_median","count":f"{gene_name}_count"}, inplace=True)
	return summary

def MAG_base_pnps(root, ocean, depths):
	bin_base = get_binning_info(root, ocean)
	bin_base = merge_single_copy_gens(bin_base, root, ocean)
	pnps = bin_base.merge(get_all_pnps(depths, root, ocean))
	pnps = merge_COG_category(pnps, root, ocean)
	whole_MAG_pnps = pnps.groupby(['bin','sample_id','depth']).agg(
		MAG_pnps_median=('pnps', np.nanmedian),
		total_count=('pnps', np.size)
	).reset_index()
	defense_pnps = pnps[pnps['category_accession'].str.contains("V", na=False)]
	defense_summary = defense_pnps[["bin","sample_id","pnps"]].groupby(['bin','sample_id']).agg(
		defense_pnps_median=('pnps', np.nanmedian),
		defense_count=('pnps', np.size)
	).reset_index()
	scg_summary = get_pnps_summary(pnps, "scg", "scg")
	out = whole_MAG_pnps.merge(scg_summary, on=["bin","sample_id"])
	# toxin_summary = get_pnps_summary(pnps, "toxin_id", "toxin")
	out = out.merge(defense_summary, on=["bin","sample_id"])
	return out

def get_res(ocean, depths):
	root="/researchdrive/zhongj2/MAG_pnps_redux"
	pnps = MAG_base_pnps(root, ocean, depths)
	ocean_sum = pnps[(pnps.defense_count >= 5) & (pnps.scg_count >= 5)] # (out.toxin_count >= 5) &
	ocean_sum = ocean_sum[(ocean_sum.MAG_pnps_median > 0) & (ocean_sum.scg_pnps_median > 0)]
	ocean_sum = ocean_sum.assign(ocean=ocean)
	for scale in ['scg','MAG']:
		ratio, denominator =f"ratio_all_{scale}", f"{scale}_pnps_median"
		ocean_sum[ratio] = ocean_sum.apply(lambda row: row.defense_pnps_median / row[denominator], axis=1)
	pnps_ratio_summary = ocean_sum[['depth', ratio]].groupby('depth').describe()[ratio][['50%','count']].reset_index()
	pnps_ratio_summary.rename(columns= {"50%":f"{ocean}_median"}, inplace=True)
	return ocean_sum, pnps_ratio_summary

def main():
	oceans = ["IN","SAT","NAT","SP","NP"]
	depths = ["SRF", "DCM", "MES", "deep"]
	per_MAG_summary, scg_defense_pnps_summary = get_res(oceans[0], depths)
	for ocean in oceans[1:]:
		pMs, os=get_res(ocean, depths)
		per_MAG_summary = per_MAG_summary.append(pMs)
		scg_defense_pnps_summary = scg_defense_pnps_summary.merge(os, on="depth")

	per_MAG_summary.to_csv(f"per_MAGs_pnps_summary.csv", index=False)
	out_f = f"per_oceans_scg_vs_defense_pnps_summary.csv"
	scg_defense_pnps_summary.to_csv(out_f, index=False)
	print(f"wrote -------- {out_f} -------- !!!")

if __name__ == "__main__":
	main()


