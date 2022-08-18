import csv
import glob
import os
import pandas as pd

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

base1="/researchdrive/zhongj2/BLAST_tara_OM-RGC_v2/"
base2="/researchdrive/zhongj2/deep_ocean_bins/deep_metagenome_transposase_BLAST"

''' a summary.signalp5 output file structure:
# SignalP-5.0   Organism: gram- Timestamp: 20220703215412
# ID    Prediction      SP(Sec/SPI)     TAT(Tat/SPI)    LIPO(Sec/SPII)  OTHER   CS Position
OM-RGC.v2.010671536     LIPO(Sec/SPII)  0.175084        0.000557        0.722455        0.101903        CS pos: 24-25. FCG-LS. Pr: 0.4868
OM-RGC.v2.010591442     OTHER   0.006607        0.001668        0.000785        0.990940
OM-RGC.v2.012952042     OTHER   0.003266        0.000160        0.000483        0.996091
'''

def get_all_signalp(signalp_file):
	out = []
	if os.path.getsize(signalp_file) == 0: # the files I created manually,
		# this MAG doesn't have any blastp match to peptidase
		return out
	signalp = list(csv.reader(open(signalp_file), delimiter='\t'))[2:]
	for line in signalp:
		Gene_ID, signal_type = line[0], line[1]
		if signal_type in ['SP(Sec/SPI)','TAT(Tat/SPI)','LIPO(Sec/SPII)']:
			out.append(Gene_ID)
	return out

# def signalp_by_domain(signalp_files):
# 	grampos_list, gramneg_list, arch_list = [], [], []
# 	for f in signalp_files:
# 		signal_list = get_all_signalp(f)
# 		if "gramPos" in f:
# 			grampos_list += signal_list
# 		elif "gramNeg" in f:
# 			gramneg_list += signal_list
# 		elif "arch" in f:
# 			arch_list += signal_list
# 		else:
# 			print(f"error. don't know what file it is: {f}")
# 			exit(2)
# 	all_signalp = set(grampos_list + gramneg_list + arch_list)
# 	domain_signalp = {"gramPos":set(grampos_list), "gramNeg":set(gramneg_list), "arch":set(arch_list)}
# 	return all_signalp, domain_signalp

def get_signalp_with_domain(signalp_files):
	out = []
	for f in signalp_files:
		signal_list = get_all_signalp(f)
		source = "error"
		if "gramPos" in f: source = "gramPos"
		elif "gramNeg" in f: source = "gramNeg"
		elif "arch" in f: source = "arch"
		else:
			print(f"error. don't know what file it is: {f}")
			exit(2)
		out+=[(gene_ID, source) for gene_ID in signal_list]
	return out

def get_OM_RGC_domain_dict(signal_CAZ_with_source, signal_pep_with_source):
	all_CAZ = set([x[0] for x in signal_CAZ_with_source])
	all_pep = set([x[0] for x in signal_pep_with_source])
	all_potential_signal = all_CAZ.union(all_pep)
	OM_RGC_f = open(f"{base1}/OM-RGC_v2.tsv", 'r')
	print("wait a while, there are 46 million lines")
	domain_dict, count = {}, 0
	while True:
		if count%3000000 == 0:
			print(f"I am at line {count}")
		count += 1
		line = OM_RGC_f.readline()
		if not line:
			break
		l = line.split("\t")
		OM_RGC_ID, domain, taxon_class = l[1], l[5], l[6]
		if OM_RGC_ID in all_potential_signal:
			if domain == "Bacteria":
				if taxon_class == "Actinobacteria" or taxon_class == "Firmicutes":
					domain_dict[OM_RGC_ID] = "gramPos"
				else:
					domain_dict[OM_RGC_ID] = "gramNeg"
			elif domain == "Archaea":
				domain_dict[OM_RGC_ID] = "arch"
			else:
				domain_dict[OM_RGC_ID] = "no match"
	OM_RGC_f.close()
	return domain_dict

''' use this to extract taxon numbers from NCBI taxonomy
tax=fullnamelineage.dmp
grep Archaea $tax | awk '{print $1 " arch"}' > TaxID_category.txt
grep -E "Actinobacteria|Firmicutes" $tax | awk '{print $1 " gramPos"}' >> TaxID_category.txt
grep Bacteria $tax | grep -vE "Actinobacteria|Firmicutes" | awk '{print $1 " gramNeg"}' >> TaxID_category.txt
'''


def get_secretory_file(outfile, signal_with_source, gene_ID_domain_dict):
	secretory_ID = []
	for tuple in signal_with_source:
		gene_ID, source = tuple
		if gene_ID_domain_dict[gene_ID] == source: # arch/gramPos/gramNeg matches!
			secretory_ID.append([gene_ID, source])
	secretory_df = pd.DataFrame(secretory_ID).astype("str")
	print(secretory_df)
	secretory_df.columns = ["gene_ID", "taxon"]
	print(secretory_df)
	secretory_df.to_csv(path_or_buf=outfile, index=False)
	if "deep" in outfile:
		secretory_df["gene_ID"].to_csv(f"{base2}/{outfile}", index=False)
	else:
		secretory_df["gene_ID"].to_csv(f"{base1}/{outfile}", index=False)

def main():
CAZ_signalp_outfiles = glob.glob("../../OM-RGC_CAZyme_signalp/*_summary.signalp5")
pep_signalp_outfiles = glob.glob("../../OM-RGC_peptidase_signalp/*_summary.signalp5")

signal_CAZ_with_source = get_signalp_with_domain(CAZ_signalp_outfiles)
signal_pep_with_source = get_signalp_with_domain(pep_signalp_outfiles)
print("start reading the big OM-RGC reference")
signal_domain_dict = get_OM_RGC_domain_dict(signal_CAZ_with_source, signal_pep_with_source)

get_secretory_file("secretory_CAZyme.txt", signal_CAZ_with_source, signal_domain_dict)
get_secretory_file("secretory_peptidase.txt", signal_pep_with_source, signal_domain_dict)
	
	# deep ocean
	deep_CAZ_signalp_f = glob.glob("../../deep_metagenome/deep_CAZyme_*_summary.signalp5")
	deep_pep_signalp_f = glob.glob("../../deep_metagenome/deep_peptidase_*_summary.signalp5")
	get_secretory_file("deep_secretory_CAZyme.txt", deep_CAZ_signalp_f)
	get_secretory_file("deep_secretory_peptidase.txt", deep_pep_signalp_f)

if __name__ == "__main__":
	main()
