#!/bin/python3
import csv
import glob
import os
from Bio import SeqIO
from first_get_CAZenzyme_and_peptidase import read_peptidase
from third_get_secretory_CAZenzyme import get_all_signalp
from MAGs_first_get_CAZenzyme_peptidase_seq import read_overview_2_agree

def clean_name(MAG_name):
	if "mp-deep_mag-" in MAG_name:
		return MAG_name.replace("mp-deep_mag-","deep_MAG_")
	else:
		return MAG_name.replace("TOBG_","")

def get_signalp_set(signalp_pos_f):
	signalp_neg_f = signalp_pos_f.replace("gramPositive", "gramNegative")
	pos_signalp_list = get_all_signalp(signalp_pos_f)
	neg_signalp_list = get_all_signalp(signalp_neg_f)
	signalp_set = set(pos_signalp_list)|set(neg_signalp_list)
	return signalp_set

def count_signalp(signalp_pos_f):
	if not os.path.exists(signalp_pos_f):
		print(f"can't find this file: {signalp_pos_f}")
	else:
		signalp_set = get_signalp_set(signalp_pos_f)
		return len(signalp_set)

def cal_length_signalp(signalp_pos_f, ORF_calls_f):
	signalp_sum_bp, ORFs_sum_bp, ORF_count = 0, 0, 0
	signalp_set = set()
	if not os.path.exists(signalp_pos_f):
		reason = "They don't have CAZyme or petidase in MAG to begin with. So there is no signalp detection"
		print(f"can't find this file: {signalp_pos_f}. {reason}")
	else:
		signalp_set = get_signalp_set(signalp_pos_f)
	fasta_seq = SeqIO.parse(open(ORF_calls_f), 'fasta')
	for f in fasta_seq:
		name, seq = f.id, str(f.seq)
		ORFs_sum_bp += len(seq)
		ORF_count += 1
		if name in signalp_set:
			signalp_sum_bp += len(seq)
	return signalp_sum_bp, ORFs_sum_bp, ORF_count

def cal_percent(d,n): return 0 if (n == 0 or not n or not d) else d/n*100

def cal_CAZenzyme(overview_f):
	signalp_pos_f = overview_f.replace("overview.txt", "two_agree_CAZenzyme_gramPositive_summary.signalp5")
	ORF_calls_f = overview_f.replace("overview.txt", "uniInput")
	signalp_count = count_signalp(signalp_pos_f)
	total_CAZ_count = len(read_overview_2_agree(overview_f))
	signalp_bp, MAG_ORFs_bp, ORF_count = cal_length_signalp(signalp_pos_f, ORF_calls_f)
	return [signalp_count, cal_percent(signalp_count,total_CAZ_count), signalp_bp, MAG_ORFs_bp, ORF_count]

def cal_peptidase(overview_f):
	signalp_pos_f = overview_f.replace("overview.txt", "diamond_peptidase_gramPositive_summary.signalp5")
	ORF_calls_f = overview_f.replace("overview.txt", "uniInput")
	signalp_count = count_signalp(signalp_pos_f)
	blastp_f = overview_f.replace("overview.txt", "peptidase_diamond_unique.blastp")
	total_pep_count = len(read_peptidase(blastp_f, split = False))
	signalp_bp, MAG_ORFs_bp, ORF_count = cal_length_signalp(signalp_pos_f, ORF_calls_f)
	return [signalp_count, cal_percent(signalp_count,total_pep_count), signalp_bp]

def get_ocean_summary_count(ocean):
	ocean_summary = []
	overviews = glob.glob(f"../../bins/{ocean}/*/overview.txt")
	for ov in overviews:
		MAG_name = clean_name(ov.split("/")[4])
		entire_row = [MAG_name] + cal_CAZenzyme(ov) + cal_peptidase(ov)
		ocean_summary.append(entire_row)
	return ocean_summary

def main():
	oceans = ["ARS","CPC","deep","EAC","IN","MED","NAT","NP","RS","SAT","SP"]
	CAZ_sum_f = "MAGs_two_agree_CAZyme_and_peptidase_signalp_summary.csv"
	colnames = ["bin",
	"signal_CAZ_count", "percent_sect_CAZ", "signal_CAZ_aa", "MAG_ORFs_aa", "ORF_count"
	"signal_pep_count", "percent_sect_pep", "signal_pep_aa"]

	with open(CAZ_sum_f, 'w') as a:
		csv_write = csv.writer(a, delimiter=',')
		csv_write.writerow(colnames)
		for ocean in oceans:
			per_ocean_summary = get_ocean_summary_count(ocean)
			csv_write.writerows(per_ocean_summary)
			print(f"summarized the -- {len(per_ocean_summary)} -- MAGs in the {ocean} Ocean")

if __name__ == "__main__":
	main()
