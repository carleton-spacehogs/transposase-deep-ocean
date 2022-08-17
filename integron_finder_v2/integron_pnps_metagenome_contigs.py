import sys
import pandas as pd
import csv
import glob
sys.path.insert(0, '../MAG_pnps_redux')
from MAG_integron_selection import explode_COG, integron_dict_helper, merge_anvi_gene_call

output_col = ["contig","ocean","gene_callers_id","is_integron","COG_accession","COG_category"]
oceans=["IN", "SAT", "NAT", "SP", "NP", "CPC", "ARS", "RS", "EAC", "MED", "deep"]

def integron_dict_of_anvio(file = "metagenome-integron-all.csv"):
	out_dict = {}
	for l in list(csv.reader(open(file)))[1:]:
		if l[-3] == "protein":
			contig, anvi_id, ocean = str(l[1]), l[0], l[-4]
			contig_key = f"{ocean}-c_{contig.split('_', 1)[1]}"
			if ocean == "deep":
				contig_key = f"{ocean}-k127_{contig.split('_', 1)[1]}"
			if contig_key in out_dict:
				out_dict[contig_key].add(anvi_id)
			else:
				out_dict[contig_key] = set(anvi_id)
	return out_dict

def read_pnps_v2(pnps_f):
	sample_pnps = pd.read_csv(pnps_f, sep = "\t").astype(str)
	sample_pnps = sample_pnps[sample_pnps['pNpS_gene_reference'].astype(float) < 100] # I just need a rank (percentile)
	sample_name = pnps_f.split("/")[-3]
	if "deep" in pnps_f:
		sample_name = (pnps_f.split("/")[-2]).split("_")[0]
	pnps_colname = "pnps_" + sample_name
	sample_pnps.rename(columns= {"corresponding_gene_call":"gene_callers_id","pNpS_gene_reference":pnps_colname}, inplace=True)
	sample_pnps = sample_pnps[["gene_callers_id", pnps_colname]]
	return sample_pnps, pnps_colname

def ocean_integron_contig(ocean):
	# ocean="IN"
	ocean="SP"
	base=f"/researchdrive/zhongj2/Tara_Oceans_Bins/{ocean}_bins_v2"
	anvio_db=f"{base}/all-gene-calls.txt"
	cog_category=f"{base}/all-COG-categories.txt"
	pnps_files=glob.glob(f"/researchdrive/zhongj2/metagenome_pnps_redux/{ocean}/profile_DBs/*/pn-ps/pNpS.txt")
	if ocean == "deep":
		base="/researchdrive/zhongj2/deep_ocean_bins"
		anvio_db=f"{base}/anvio-gene-call.txt"
		cog_category=f"{base}/all-anvi-cog-category.tsv"
		pnps_files=glob.glob(f"{base}/per-sample-mapping/*_pn-ps/pNpS.txt")
	integron_dict = integron_dict_of_anvio()
	all_genes = list(csv.reader(open(anvio_db), delimiter='\t'))[1:]
	integron_contig = []
	for l in all_genes:
		anvi_id, contig = l[0], l[1]
		contig_key = f"{ocean}-{contig}"
		if contig_key in integron_dict:
			if anvi_id in integron_dict[contig_key]:
				integron_contig.append([contig, ocean, anvi_id, "integron"])
			else:
				integron_contig.append([contig, ocean, anvi_id, "non-integron"])
	integron_contig_df = pd.DataFrame(integron_contig)
	integron_contig_df.columns = output_col[:4]
	COG_category = pd.read_csv(cog_category, sep = "\t").astype(str)
	integron_contig_df = integron_contig_df.merge(COG_category[['gene_callers_id',
		'accession','function']], on='gene_callers_id', how = "left")
	integron_contig_df.columns = output_col
	samples = []
	for pnps_f in pnps_files:
		pnps, pnps_colname = read_pnps_v2(pnps_f)
		samples.append(pnps_colname)
		integron_contig_df = integron_contig_df.merge(pnps, on='gene_callers_id', how = "left")
	integron_contig_df = integron_contig_df.dropna(subset=samples, how='all')
	integron_contig_df.to_csv(f'per_contig_integron/{ocean}.csv', index=False)
	print(f"I got a total of: {len(set(integron_contig_df.contig.values.tolist()))} integrons in {ocean}")

# ocean_integron_contig("IN")
# ocean_integron_contig("deep")
ocean_integron_contig("SP")

