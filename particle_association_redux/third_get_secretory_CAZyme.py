import csv
import glob
import os
import pandas as pd

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

base1="/researchdrive/zhongj2/BLAST_tara_OM-RGC_v2"
base2="/researchdrive/zhongj2/deep_ocean_bins"

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
	OM_RGC_f.close()
	return domain_dict

''' use run_dbcan to convert nucleotide into proteins
base="/researchdrive/zhongj2/deep_ocean_bins"
db=${base}/OM-RGC-reference/geneDB
query=${base}/deep_metagenome_transposase_BLAST/deep_anvi-genecall.fasta
# sed -i 's/_1 / /g' ${db}.faa
# diamond makedb --in ${db}.faa -d $db -p 10
diamond blastx -d $db -q $query --threads 20 --evalue 0.00001 --more-sensitive -o ${base}/OM-RGC-reference/db_vs_anvio.blastx
cat db_vs_anvio.blastx | sort -k11 | sort -u -k1,1 > db_vs_anvio_unique2.blastx
'''

def deep_OM_RGC_domain_base():
	db_file = f"{base2}/OM-RGC-reference/M-geneDB.v1_1115269-seqs_full-header-annotation.txt"
	db = list(csv.reader(open(db_file), delimiter='_'))
	out = []
	for l in db:
		gene_ID = f"{l[0]}_{l[1]}"
		if len(l) > 11:
			if "Archaea" in l[10]: 
				out.append([gene_ID, "arch"])
			elif "Bacteria" in l[10]:
				if l[11] == "Actinobacteria" or l[11] == "Firmicutes":
					out.append([gene_ID, "gramPos"])
				else:
					out.append([gene_ID, "gramNeg"])
	ref_source = pd.DataFrame(out)
	ref_source.columns = ["ref","source"]
	return ref_source

def get_deep_OM_RGC_domain(): # now deprecated
	ref_source = deep_OM_RGC_domain_base()
	domain_dict = {}
	blast_match = pd.read_csv(f"{base2}/OM-RGC-reference/db_vs_anvio_unique2.blastx", sep="\t", header=None)
	blast_match = blast_match[[0,1]]
	blast_match.columns = ["anvi_id", "ref"]
	ref_source = ref_source.merge(blast_match, on = "ref")
	for l in ref_source[["anvi_id","source"]].values.tolist():
		domain_dict[str(l[0])] = l[1]
	return domain_dict


def get_secretory_file(outfile, signal_with_source, gene_ID_domain_dict):
	secretory_ID = []
	for tuple in signal_with_source:
		gene_ID, source = tuple
		if gene_ID in gene_ID_domain_dict and gene_ID_domain_dict[gene_ID] == source: # arch/gramPos/gramNeg matches!
			secretory_ID.append([gene_ID, source])
	secretory_df = pd.DataFrame(secretory_ID).astype("str")
	secretory_df.columns = ["gene_ID", "taxon"]
	secretory_df.to_csv(path_or_buf=outfile, index=False)
	if "deep" in outfile:
		secretory_df["gene_ID"].to_csv(f"{base2}/deep_metagenome_transposase_BLAST/{outfile}", index=False)
	else:
		secretory_df["gene_ID"].to_csv(f"{base1}/{outfile}", index=False)

def filter_deep_CAZpep_without_taxon(old_CAZpep, norm_CAZpep, domain_dict):
	gene_id = set(line.strip() for line in open(old_CAZpep))
	print(f"done reading: {old_CAZpep} ; total: {len(gene_id)} lines")
	out_lines, count = "", 0
	for id in gene_id:
		count += 1
		if count % 1000 != 0:
			print(f"I am at line {count}")
		if id in domain_dict:
			out_lines += (str(gene_id) + "\n")
	with open(norm_CAZpep, 'w') as out_file:
		out_file.write(out_lines)

# # tooooooooo slow
# def filter_deep_CAZpep_without_taxon(old_CAZpep, norm_CAZpep, domain_dict):
# 	gene_id = list(line.strip() for line in open(old_CAZpep))
# 	gene_id_df = pd.DataFrame(gene_id)
# 	print(f"done reading: {old_CAZpep} ; total: {len(gene_id)} lines")
# 	out_lines, count = "", 0
# 	for id in gene_id:
# 		count += 1
# 		if count % 1000 != 0:
# 			print(f"I am at line {count}")
# 		if id in domain_dict:
# 			out_lines += (str(gene_id) + "\n")
# 	with open(norm_CAZpep, 'w') as out_file:
# 		out_file.write(out_lines)

def filter_deep_CAZpep_without_taxon(old_CAZpep, norm_CAZpep, domain_dict):
	gene_id_df = pd.read_csv(old_CAZpep, header=None).astype(str)
	avail_df = pd.DataFrame(list(domain_dict.keys())).astype(str)
	print("done reading files, merging")
	intersect = gene_id_df.merge(avail_df, on=0, how="inner")
	intersect.to_csv(norm_CAZpep, index=False, header=None)

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
	# deep_CAZ_signalp_f = glob.glob("../../deep_metagenome/deep_CAZyme_*_summary.signalp5")
	# deep_pep_signalp_f = glob.glob("../../deep_metagenome/deep_peptidase_*_summary.signalp5")
	# deep_signal_CAZ = get_signalp_with_domain(deep_CAZ_signalp_f)
	# deep_signal_pep = get_signalp_with_domain(deep_pep_signalp_f)
	# deep_domain_dict = get_deep_OM_RGC_domain() # it's a small reference, don't need to select

	deep_CAZ_signalp_f = glob.glob("../../deep-OM-RGC/CAZyme_*_summary.signalp5")
	deep_pep_signalp_f = glob.glob("../../deep-OM-RGC/peptidase_*_summary.signalp5")
	deep_signal_CAZ = get_signalp_with_domain(deep_CAZ_signalp_f)
	deep_signal_pep = get_signalp_with_domain(deep_pep_signalp_f)
	deep_domain = deep_OM_RGC_domain_base()
	deep_dom_dict = { l[0]:l[1] for l in deep_domain[["ref","source"]].values.tolist() }

	get_secretory_file("deep_secretory_CAZyme.txt", deep_signal_CAZ, deep_dom_dict)
	get_secretory_file("deep_secretory_peptidase.txt", deep_signal_pep, deep_dom_dict)

	'''some secretory CAZyme does not have a taxon origin, so they were not counted as signal molecule
	need to account for the lost that from the total CAZyme, too. '''
	# already done under the new OM-RGC based reference
	# base3=f"{base2}/deep_metagenome_transposase_BLAST"
	# filter_deep_CAZpep_without_taxon(f"{base3}/CAZyme/deep_CAZyme_anvigeneID.txt", f"{base3}/all-CAZyme-norm.txt", deep_domain_dict)
	# filter_deep_CAZpep_without_taxon(f"{base3}/peptidase/deep_peptidase_anvigeneID.txt", f"{base3}/all-peptidase-norm.txt", deep_domain_dict)

	filter_deep_CAZpep_without_taxon(f"{base1}/all_CAZyme.txt", f"{base1}/all_CAZyme_norm.txt", signal_domain_dict)
	filter_deep_CAZpep_without_taxon(f"{base1}/all_peptidase.txt", f"{base1}/all_peptidase_norm.txt", signal_domain_dict)

if __name__ == "__main__":
	main()
