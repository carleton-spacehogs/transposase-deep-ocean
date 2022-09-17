import cal_cov_tara as tara
import csv
import pandas as pd

base="/researchdrive/zhongj2/deep_ocean_bins"
base1=f"{base}/OM-RGC-ref"

def read_blastp(f, col):
	res = pd.read_csv(f, sep = "\t")
	gene_IDs = pd.DataFrame(res.iloc[:,col])
	gene_IDs.columns = ["gene_ID"]
	return gene_IDs

def get_gene_specific_ref_through_COG(search_string):
	cogfile = pd.read_csv('../COG_numbers_categories_annotations.txt', sep = "\t")
	cogf = cogfile.iloc[:,[0,2]]
	cogf.columns = ['COG_num', "COG_category"]
	COG_nums = cogf[cogf.COG_category.str.contains(search_string)].COG_num.values.tolist()
	COG_cat_nums = "|".join(COG_nums)
	return COG_cat_nums

def cal_sum_merge_md(gene_rows, md, gene):
	gene = gene.replace("transport_and_metabolism", "TM")
	gene = f"DNA_{gene}"
	print(gene)
	gene_abun = gene_rows.filter(like='SRR', axis=1)
	sample_abun = pd.DataFrame(gene_abun.sum()).reset_index()
	sample_abun.columns = ["sample", gene]
	md = md.merge(sample_abun, on = "sample")
	md[gene] = md[gene]*150/1000000000/md.Sequencing_depth_Gbp
	return md

def main():
	reference = pd.read_csv(f"{base1}/deep-ref-abundance.csv")
	reference = reference.rename(columns = {'id':'gene_ID'})

	biofilm_id = read_blastp(f"{base1}/diamond-blastp/biofilm_diamond.blastp", 1)
	trans_id = read_blastp(f"{base1}/diamond-blastp/transposase_diamond.blastp", 1)
	CAZ_id = pd.read_csv("../particle_association_redux/deep_CAZyme_ID.txt", names=["gene_ID"], header=None)
	pep_id = pd.read_csv("../particle_association_redux/deep_peptidase_ID.txt", names=["gene_ID"], header=None)
	sect_CAZ_id = pd.read_csv("../particle_association_redux/deep_secretory_CAZyme.txt")
	sect_pep_id = pd.read_csv("../particle_association_redux/deep_secretory_peptidase.txt")

	md_cols = "sample,station,lower_filter_size,Ocean_DNA,Ocean,Depth,Longitude,Latitude,Salinity_PSU,Temperature,Oxygen_DNA,WATERMASS,Basin,Sequencing_depth_Gbp".split(",")
	metadata = pd.read_csv(f"{base}/deep_metagenome_transposase_BLAST/sample_metadata.tsv",
						sep = "\t", names = md_cols, header = 0)

	single_genes_id = [biofilm_id, trans_id, CAZ_id, pep_id, sect_CAZ_id, sect_pep_id]
	for i in range(len(tara.single_genes)):
		gene_rows = single_genes_id[i].merge(reference, on="gene_ID")
		metadata = cal_sum_merge_md(gene_rows, metadata, tara.single_genes[i])

	for cog_category in tara.COG_list:
		COG_cat_id = get_gene_specific_ref_through_COG(cog_category)
		gene_rows = reference[reference.COG.str.contains(COG_cat_id)]
		metadata = cal_sum_merge_md(gene_rows, metadata, cog_category)

	metadata.to_csv("Malaspina-genes-coverage-OM-RGC.csv", index=False)

if __name__ == "__main__":
	main()