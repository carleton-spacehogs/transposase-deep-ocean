import csv
import glob
from third_get_secretory_CAZenzyme import get_all_signalp
from MAGs_first_get_CAZenzyme_seq import read_overview

def MAGs_sectory_gene_count(signalp_pos_f, gene):
	# signalp_pos_f, e.g.: ../../bins/ARS/TOBG_ARS-1426/CAZenzyme_gramPositive_summary.signalp5
	signalp_neg_f = signalp_pos_f.replace("gramPositive", "gramNegative")
	overview_f = signalp_pos_f.replace(f"{gene}_gramPositive_summary.signalp5", "overview.txt")
	pos_CAZ_list = get_all_signalp(signalp_pos_f)
	neg_CAZ_list = get_all_signalp(signalp_neg_f)
	signal_CAZ_set = set(pos_CAZ_list)|set(neg_CAZ_list)
	signal_CAZ_count = len(signal_CAZ_set)
	total_CAZ_count = len(read_overview(overview_f))
	MAG_name = signalp_pos_f.split("/")[4]
	return [MAG_name, signal_CAZ_count, total_CAZ_count]

def get_ocean_summary_count(ocean, gene):
	ocean_summary = []
	ocean_gene_signalp_pos = glob.glob(f"../../bins/{ocean}/*/{gene}_gramPositive_summary.signalp5")
	for signalp_pos_f in ocean_gene_signalp_pos:
		MAG_summary = MAGs_sectory_gene_count(signalp_pos_f, gene)
		ocean_summary.append(MAG_summary)
	return ocean_summary

def main():
	# oceans = ["ARS","CPC","deep","EAC","IN","MED","NAT","NP","RS","SAT","SP"]
	oceans = ["ARS", "CPC", "IN", "NP", "RS"]
	CAZ_sum_f = "MAGs_CAZenzyme_signalp_summary.csv"
	sum_f_col = "bin,signal_CAZ_count,total_CAZ_count\n"
	with open(CAZ_sum_f, 'w') as a:
		a.write(sum_f_col)
		csv_write = csv.writer(a, delimiter=',')
		for ocean in oceans:
			per_ocean_summary = get_ocean_summary_count(ocean, "CAZenzyme")
			csv_write.writerows(per_ocean_summary)
			print(f"summarized the -- {len(per_ocean_summary)} -- MAGs in the {ocean} Ocean")

if __name__ == "__main__":
	main()
