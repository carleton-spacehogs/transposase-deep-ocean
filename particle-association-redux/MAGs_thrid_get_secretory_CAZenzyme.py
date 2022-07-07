import csv
import glob
from third_get_secretory_CAZenzyme import get_all_signalp
from MAGs_first_get_CAZenzyme_seq import read_overview

def MAGs_sectory_gene_count(signalp_pos_f, signalp_output_file, include_MAG_name = True):
	# signalp_pos_f, e.g.: ../../bins/ARS/TOBG_ARS-1426/diamond_CAZenzyme_gramPositive_summary.signalp5
	signalp_neg_f = signalp_pos_f.replace("gramPositive", "gramNegative")
	overview_f = signalp_pos_f.replace(signalp_output_file, "overview.txt")
	pos_signalp_list = get_all_signalp(signalp_pos_f)
	neg_signalp_list = get_all_signalp(signalp_neg_f)
	signalp_set = set(pos_signalp_list)|set(neg_signalp_list)
	signalp_count = len(signalp_set)
	total_signalp_count = len(read_overview(overview_f))
	if include_MAG_name:
		return [signalp_pos_f.split("/")[4], signalp_count, total_signalp_count]
	else:
		return [signalp_count, total_signalp_count]

def get_ocean_summary_count(ocean, CAZenzyme_signalp, peptidase_signalp):
	ocean_summary = []
	CAZenzyme_signalp_GramPos = glob.glob(f"../../bins/{ocean}/*/{CAZenzyme_signalp}")
	peptidase_signalp_GramPos = glob.glob(f"../../bins/{ocean}/*/{peptidase_signalp}")
	if len(CAZenzyme_signalp_GramPos) == len(peptidase_signalp_GramPos):
		for i in len(CAZenzyme_signalp_GramPos):
			left = MAGs_sectory_gene_count(CAZenzyme_signalp_GramPos[i], CAZenzyme_signalp) # I can get away with only 1 argument, the second argu just make my code easier to write
			right =  MAGs_sectory_gene_count(peptidase_signalp_GramPos[i], peptidase_signalp, include_MAG_name = False)
			entire_row = left + right
			ocean_summary.append(entire_row)
	else:
		exit("Each MAGs show have a gramPositive_summary.signalp5 for both CAZenzyme and peptidase")
	return ocean_summary

def main():
	# oceans = ["ARS","CPC","deep","EAC","IN","MED","NAT","NP","RS","SAT","SP"]
	oceans = ["ARS", "CPC", "IN", "NP", "RS"]
	CAZ_sum_f = "MAGs_CAZenzyme_and_peptidase_signalp_summary.csv"
	sum_f_col = "bin,signal_CAZ_count,total_CAZ_count\n"
	with open(CAZ_sum_f, 'w') as a:
		a.write(sum_f_col)
		csv_write = csv.writer(a, delimiter=',')
		for ocean in oceans:
			per_ocean_summary = get_ocean_summary_count(ocean,
			"diamond_CAZenzyme_gramPositive_summary.signalp5",
			"diamond_peptidase_gramPositive_summary.signalp5")
			csv_write.writerows(per_ocean_summary)
			print(f"summarized the -- {len(per_ocean_summary)} -- MAGs in the {ocean} Ocean")

if __name__ == "__main__":
	main()
