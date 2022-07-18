import csv
import pandas as pd
from MAG_singleCopyGene_defense_pnps import MAG_base_pnps, get_binning_info
import numpy as np
from Bio import SeqIO

def merge_transposase(df, root, ocean):
	trans=pd.read_csv(f"{root}/{ocean}/transposase_diamond_unique.blastp", sep = "\t", header=None).astype(str)[[1,0]]
	trans.rename(columns={1:"gene_callers_id", 0:"transposase_id"}, inplace=True)
	df = df.merge(trans, on="gene_callers_id", how = "left")
	return df

def get_gene_length(root, ocean):
	in_faa = f"{root}/{ocean}/all_bins_db/anvi-gene-calls.faa"
	fasta_seq = SeqIO.parse(open(in_faa), 'fasta')
	out = []
	for f in fasta_seq:
		gene_caller_id, seq_len = f.id, len(str(f.seq))
		out.append([gene_caller_id, seq_len])
	gene_length = pd.DataFrame(out)
	gene_length.columns = ["gene_callers_id", "length"]
	return gene_length

def gene_abun_base(ocean, depth, root):
	abun_file = f"{root}/{ocean}/PROFILE-{depth}/all-GENE-COVERAGES.txt"
	abun = pd.read_csv(abun_file, sep = "\t")
	abun.rename(columns={"key":"gene_callers_id"}, inplace=True)
	tmp = abun.drop(['gene_callers_id'], axis=1)
	abun = abun[~(tmp < 0.01).all(axis=1)].astype(str)
	abun = abun.merge(get_binning_info(root, ocean)[["gene_callers_id","bin"]], on = "gene_callers_id", how = "left")
	abun = abun.merge(get_gene_length(root, ocean), on = "gene_callers_id", how = "left")
	return abun

def cal_MAG_abundance(ocean, depth, root):
	abun = gene_abun_base(ocean, depth, root)
	abun_bp = abun.filter(regex='ERR|SRR').astype(float).multiply(abun.length, axis="index")
	abun_bp2 = abun[["bin"]].join(abun_bp)
	bin_sum = abun_bp2.groupby(['bin']).agg('sum').reset_index()
	f_name = f"MAGs-coverage-{ocean}-{depth}.csv"
	bin_sum.to_csv(f_name, index=False)
	print(f"!!! -------- wrote {f_name} ----------")

def cal_transposase_per_MAG(ocean, depth, root):
	trans=pd.read_csv(f"{root}/{ocean}/transposase_diamond_unique.blastp", sep = "\t", header=None).astype(str)[[1,0]]
	trans.rename(columns={1:"gene_callers_id", 0:"transposase_id"}, inplace=True)
	abun = gene_abun_base(ocean, depth, root)
	trans = trans.merge(abun, on = "gene_callers_id")
	trans_bp = trans.filter(regex='ERR|SRR').astype(float).multiply(trans.length, axis="index")
	trans_bp2 = trans[["bin"]].join(trans_bp)
	trans_sum = trans_bp2.groupby(['bin']).agg('sum').reset_index()
	MAG_cov = pd.read_csv(f"MAGs-coverage-{ocean}-{depth}.csv")
	MAG_cov = trans_sum[["bin"]].merge(MAG_cov, on = "bin").reset_index()
	prop_trans = trans_sum.filter(regex='ERR|SRR')/MAG_cov.filter(regex='ERR|SRR')
	prop_trans = trans_sum[['bin']].join(prop_trans)
	f_name = f"prop-transposase-{ocean}-{depth}.csv"
	prop_trans.to_csv(f_name, index=False)
	print(f"!!! -------- wrote {f_name} ----------")

def trans_scg_per_ocean(root, ocean, depths):
	prop_trans = pd.read_csv(f"prop-transposase-{ocean}-{depths[0]}.csv")
	for depth in depths[1:]:
		pt = pd.read_csv(f"prop-transposase-{ocean}-{depth}.csv")
		prop_trans = prop_trans.merge(pt, on = "bin", how = "outer")
	prop_trans = prop_trans.fillna(0)
	pnps = MAG_base_pnps(root, ocean, depths)
	pnps = pnps[pnps.scg_count >= 5]
	transposase_abundance_col = []
	for index, row in pnps.iterrows():
		MAG_abun = prop_trans[prop_trans.bin == row.bin][row.sample_id]
		if MAG_abun.empty:
			transposase_abundance_col.append(float(-1))
		else:
			transposase_abundance_col.append(float(MAG_abun))
	pnps = pnps.assign(Trans_abun=transposase_abundance_col)
	return pnps

def main():
	oceans = ["IN","SAT","NAT","SP","NP"]
	depths = ["SRF", "DCM", "MES", "deep"]
	root="/researchdrive/zhongj2/MAG_pnps_redux"
	for ocean in oceans:
		for depth in depths:
			# cal_MAG_abundance(ocean, depth, root)
			cal_transposase_per_MAG(ocean, depth, root)

	base = trans_scg_per_ocean(root, oceans[0], depths)
	for ocean in oceans[1:]:
		base = base.append(trans_scg_per_ocean(root, ocean, depths))

	f_name = f"transposase-abun-vs-scg-pnps.csv"
	base.to_csv(f_name, index=False)
	print(f"!!!--------- wrote {f_name} ----------")


if __name__ == "__main__":
	main()



