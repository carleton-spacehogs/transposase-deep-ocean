import sys
import pandas as pd
sys.path.insert(0, '../particle_association_redux')
from MAGs_first_get_CAZenzyme_peptidase_seq import read_overview_2_agree, read_overview_only_diamond
from MAGs_third_get_secretory_summary import get_signalp_set
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, MAG_base_pnps
from MAG_singleCopyGene_vs_transposaseAbundance import gene_abun_base

# change read_overview_2_agree -> read_overview_only_diamond
# prop-CAZyme-{ocean}.csv -> prop-diamond-only-CAZyme-{ocean}.csv
# CAZyme-abun-vs-defense-pnps.csv -> diamond-only-CAZyme-abun-vs-defense-pnps.csv

# this is useless!!!! I am not using CAZyme directly
def merge_CAZyme(df, ocean):
	CAZyme_overview = f"/workspace/data/zhongj/MAG_CAZymes/{ocean}_CAZyme/overview.txt"
	gene_id_set = read_overview_only_diamond(CAZyme_overview)
	mask = [gene_id in gene_id_set for gene_id in df.gene_callers_id]
	return df[mask]

def merge_signalp_CAZyme(df, ocean):
	signal_CAZyme_pos = f"/workspace/data/zhongj/MAG_CAZymes/{ocean}-CAZenzyme_gramNeg_summary.signalp5"
	signal_CAzyme_set = get_signalp_set(signal_CAZyme_pos)
	mask = [gene_id in signal_CAzyme_set for gene_id in df.gene_callers_id]
	return df[mask]

def cal_CAZyme_per_MAG(ocean, depth, root):
	all_gene_abun = gene_abun_base(ocean, depth, root, binning_delim = " ")
	signal_CAZyme_abun = merge_signalp_CAZyme(all_gene_abun, ocean)
	signal_CAZyme_abun.groupby(['bin']).agg('sum').reset_index()
	CAZyme_bp = signal_CAZyme_abun.filter(regex='ERR|SRR').astype(float).multiply(signal_CAZyme_abun.length, axis="index")
	CAZyme_bp2 = signal_CAZyme_abun[["bin"]].join(CAZyme_bp)
	CAZyme_sum = CAZyme_bp2.groupby(['bin']).agg('sum').reset_index().reset_index()
	MAG_cov = pd.read_csv(f"intermediate-files/MAGs-coverage-{ocean}-{depth}.csv")
	MAG_cov = CAZyme_sum[["bin"]].merge(MAG_cov, on = "bin").reset_index()
	prop_CAZyme = CAZyme_sum.filter(regex='ERR|SRR')/MAG_cov.filter(regex='ERR|SRR')
	prop_CAZyme = CAZyme_sum[['bin']].join(prop_CAZyme)
	return prop_CAZyme

def loop_depth_CAZyme(ocean, depths, root):
	out = cal_CAZyme_per_MAG(ocean, depths[0], root)
	for depth in depths[1:]:
		out = out.merge(cal_CAZyme_per_MAG(ocean, depth, root), on="bin", how="outer")
	out = out.fillna(0)
	f_name = f"intermediate-files/prop-signal-CAZyme-{ocean}.csv"
	out.to_csv(f_name, index=False)
	print(f"!!! -------- wrote {f_name} ----------")

def defense_pnps_vs_CAZyme_per_ocean(root, ocean, depths):
	prop_signal_CAZyme = pd.read_csv(f"intermediate-files/prop-signal-CAZyme-{ocean}.csv")
	pnps = MAG_base_pnps(root, ocean, depths)
	# pnps = pnps[pnps.defense_count >= 5]
	signal_CAZyme_abundance_col = []
	for index, row in pnps.iterrows():
		MAG_abun = prop_signal_CAZyme[prop_signal_CAZyme.bin == row.bin][row.sample_id]
		if MAG_abun.empty:
			signal_CAZyme_abundance_col.append(float(-1))
		else:
			signal_CAZyme_abundance_col.append(float(MAG_abun))
	pnps = pnps.assign(signal_CAZyme_abun=signal_CAZyme_abundance_col)
	return pnps

def main():
	o_depths = ocean_depths()
	root="/researchdrive/zhongj2/MAG_pnps_redux"
	for ocean, depths in o_depths.items():
		print(ocean)
		loop_depth_CAZyme(ocean, depths, root)
	
	base=pd.DataFrame()
	for ocean, depths in o_depths.items():
		print(ocean)
		base = base.append(defense_pnps_vs_CAZyme_per_ocean(root, ocean, depths))
	f_name = f"singalCAZyme-abun-vs-genes-pnps.csv"
	base.to_csv(f_name, index=False)
	print(f"!!!--------- wrote {f_name} ----------")

if __name__ == "__main__":
	main()

