#!/bin/python3
import csv
import glob
import os
from Bio import SeqIO
from first_get_CAZenzyme_and_peptidase import read_peptidase
from third_get_secretory_CAZenzyme import get_all_signalp
from MAGs_first_get_CAZenzyme_peptidase_seq import read_overview_2_agree, read_overview_only_diamond

CAZyme_signal_pos = "two_agree_CAZenzyme_gramPositive_summary.signalp5"
peptidase_signal_pos = "diamond_peptidase_gramPositive_summary.signalp5"

def find_gram_positive_bact():
	'''Actinobacteria- and Firmicutes-affiliated sequences were predicted under Gram-positive mode, while other bacterial sequences were predicted under Gram-negative mode. For archaeal sequences, PSORTb (3.0.2) (57) was used to predict the subcellular location, because SignalP (4.0) (56) does not support sequences of Archaea. (Zhao et al, 2020)'''
	MAG_phylum = {}
	base = '../Figure Generating/data'
	tara = list(csv.reader(open(f'{base}/bin_taxon.csv'), delimiter=','))[1:]
	for l in tara:
		MAG_name, domain, phylum = l[0], l[15], l[16]
		MAG_phylum[MAG_name] = [domain, phylum]
	deep = list(csv.reader(open(f'{base}/malaspina_bin_taxon_over70complete.csv'), delimiter=','))[1:]
	for l in deep:
		MAG_name, domain, phylum = l[0].replace("mp-deep_mag-","deep_MAG_"), l[1], l[2]
		MAG_phylum[MAG_name] = [domain, phylum]
	return MAG_phylum

def clean_name(MAG_name):
	if "mp-deep_mag-" in MAG_name:
		return MAG_name.replace("mp-deep_mag-","deep_MAG_")
	else:
		return MAG_name.replace("TOBG_","")

def give_all(signalp_pos_f):
	signalp_neg_f = signalp_pos_f.replace("gramPositive", "gramPositive")
	pos_signalp_list = get_all_signalp(signalp_pos_f)
	neg_signalp_list = get_all_signalp(signalp_neg_f)
	signalp_set = set(pos_signalp_list)|set(neg_signalp_list)
	return signalp_set

def get_signalp_set(signalp_pos_f, MAG_name, MAG_phylum_dict):
	domain, phylum = MAG_phylum_dict[MAG_name]
	if domain == "Bacteria":
		if phylum == "Actinobacteria" or phylum == "Firmicutes":
			return set(get_all_signalp(signalp_pos_f))
		else:
			signalp_neg_f = signalp_pos_f.replace("gramPositive", "gramPositive")
			return set(get_all_signalp(signalp_neg_f))
	elif domain == "Archaea":
		pos_f = signalp_pos_f.split('/')[-1]
		type_f = "peptidase" if "peptidase" in pos_f else "CAZyme"
		psort_f = signalp_pos_f.replace(pos_f, f'signal-{type_f}_psortb_archaea.txt')
		pred = list(csv.reader(open(psort_f), delimiter='\t'))[1:]
		signal_CAZyme=set()
		for l in pred:
			if "Extracellular" == l[1]:
				signal_CAZyme.add(l[0].strip())
		print(f"Archaea MAG: {MAG_name} : {len(signal_CAZyme)} extracellular CAZyme")
		return signal_CAZyme
	else:
		print(f"Non bacteria or archaea, it was: {domain}, just give you any predict I have")
		return give_all(signalp_pos_f)

def count_signalp(signalp_pos_f, MAG_name, MAG_phylum_dict):
	if not os.path.exists(signalp_pos_f):
		print(f"can't find this file: {signalp_pos_f}")
		return 0
	else:
		signalp_set = get_signalp_set(signalp_pos_f, MAG_name, MAG_phylum_dict)
		return len(signalp_set)

def cal_length_signalp(signalp_pos_f, ORF_calls_f, MAG_name, MAG_phylum_dict):
	signalp_sum_bp, ORFs_sum_bp, ORF_count = 0, 0, 0
	signalp_set = set()
	if not os.path.exists(signalp_pos_f):
		reason = "They don't have CAZyme or petidase in MAG to begin with. So there is no signalp detection"
		print(f"can't find this file: {signalp_pos_f}. {reason}")
	else:
		signalp_set = get_signalp_set(signalp_pos_f, MAG_name, MAG_phylum_dict)
	fasta_seq = SeqIO.parse(open(ORF_calls_f), 'fasta')
	for f in fasta_seq:
		name, seq = f.id, str(f.seq)
		ORFs_sum_bp += len(seq)
		ORF_count += 1
		if name in signalp_set:
			signalp_sum_bp += len(seq)
	return signalp_sum_bp, ORFs_sum_bp, ORF_count

def cal_percent(d,n): return 0 if (n == 0 or not n or not d) else d/n*100

def cal_CAZyme(overview_f, MAG_name, MAG_phylum_dict):
	signalp_pos_f = overview_f.replace("overview.txt", CAZyme_signal_pos)
	ORF_calls_f = overview_f.replace("overview.txt", "uniInput")
	signalp_count = count_signalp(signalp_pos_f, MAG_name, MAG_phylum_dict)
	# total_CAZ_count = len(read_overview_2_agree(overview_f))
	total_CAZ_count = len(read_overview_only_diamond(overview_f))
	signalp_bp, MAG_ORFs_aa, ORF_count = cal_length_signalp(signalp_pos_f, ORF_calls_f, MAG_name, MAG_phylum_dict)
	return [signalp_count, cal_percent(signalp_count,total_CAZ_count), signalp_bp, MAG_ORFs_aa, ORF_count]

def cal_peptidase(overview_f, MAG_name, MAG_phylum_dict):
	signalp_pos_f = overview_f.replace("overview.txt", peptidase_signal_pos)
	ORF_calls_f = overview_f.replace("overview.txt", "uniInput")
	signalp_count = count_signalp(signalp_pos_f, MAG_name, MAG_phylum_dict)
	blastp_f = overview_f.replace("overview.txt", "peptidase_diamond_unique.blastp")
	total_pep_count = len(read_peptidase(blastp_f, split = False))
	signalp_bp, MAG_ORFs_bp, ORF_count = cal_length_signalp(signalp_pos_f, ORF_calls_f, MAG_name, MAG_phylum_dict)
	return [signalp_count, cal_percent(signalp_count,total_pep_count), signalp_bp]

def get_ocean_summary_count(ocean, MAG_phylum_dict):
	ocean_summary = []
	overviews = glob.glob(f"../../bins/{ocean}/*/overview.txt")
	# archaea_file = open('archaea-MAGs.txt', 'a')
	for ov in overviews:
		MAG_name = clean_name(ov.split("/")[4])
		if MAG_name in MAG_phylum_dict:
			domain = MAG_phylum_dict[MAG_name][0]
			# 	archaea_file.write(MAG_name + '\n')
			entire_row = [MAG_name, domain] + cal_CAZyme(ov, MAG_name, MAG_phylum_dict) + cal_peptidase(ov, MAG_name, MAG_phylum_dict)
			ocean_summary.append(entire_row)
		else:
			print(f"MAG {MAG_name} doesn't meet the >70 completeness and <10 containmination criteria")
	# archaea_file.close()
	return ocean_summary

def main():
	oceans = ["ARS","CPC","deep","EAC","IN","MED","NAT","NP","RS","SAT","SP"]
	CAZ_sum_f = "MAGs_diamond_CAZyme_and_peptidase_signalp_summary.csv"
	colnames = ["bin", "domain",
	"signal_CAZ_count", "percent_sect_CAZ", "signal_CAZ_aa", "MAG_ORFs_aa", "ORF_count",
	"signal_pep_count", "percent_sect_pep", "signal_pep_aa"]
	MAG_phylum_dict = find_gram_positive_bact()

	with open(CAZ_sum_f, 'w') as a:
		csv_write = csv.writer(a, delimiter=',')
		csv_write.writerow(colnames)
		for ocean in oceans:
			per_ocean_summary = get_ocean_summary_count(ocean, MAG_phylum_dict)
			csv_write.writerows(per_ocean_summary)
			print(f"summarized the -- {len(per_ocean_summary)} -- MAGs in the {ocean} Ocean")

if __name__ == "__main__":
	main()
