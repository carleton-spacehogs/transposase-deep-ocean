import pandas as pd
import utils

base = "/researchdrive/zhongj2/MAG_pnps_redux"

class bin_sample_ranker:
	def __init__(self, sample_id, bin, ocean, depth):
		self.sample_id = sample_id
		self.bin = bin
		self.ocean = ocean
		self.depth = depth
		self.init_pnps_ranking()
	def init_pnps_ranking(self):
		# find the list of pnps
		pnps_file = f"{base}/{self.ocean}/PROFILE-{self.depth}/all_MAGs_pnps/pNpS.txt"
		MAG_origin = f"{base}/{self.ocean}/all_bins_db/gene_callers_id-contig.txt"
		sample_pnps_id = self.sample_id
		if "deep" in self.bin:
			sample_pnps_id = f"{self.ocean}_{self.sample_id}"
			pnps_file = f"{base}/deep/PROFILE-{self.ocean}/all_MAGs_pnps/pNpS.txt"
			MAG_origin = f"{base}/deep/all_bins_db/gene_callers_id-contig.txt"
		pnps_file = pd.read_csv(pnps_file, sep = "\t")
		pnps_file.rename(columns= {"corresponding_gene_call":"gene_callers_id"}, inplace=True)
		MAG_origin = utils.get_bin_helper(MAG_origin)
		bin_anvi_id = MAG_origin[MAG_origin.bin == self.bin]
		sample_pnps = pnps_file[pnps_file.sample_id == sample_pnps_id].astype(str)
		bin_sample_pnps = sample_pnps.merge(bin_anvi_id, how="inner", on="gene_callers_id")
		pnps_list = bin_sample_pnps.pNpS_gene_reference.values.tolist()
		self.ranker = utils.ranking(pnps_list)
	def find_pnps_ranking(self, sample_id, bin, this_pnps):
		if sample_id == self.sample_id and bin == self.bin:
			return self.ranker.find_pnps_ranking(this_pnps)

integrons = pd.read_csv("MAG_integron_pnps_all_v2.csv")

ORF_info = ["gene_callers_id","pnps","COG_accession","COG_category"]
bin_info = ["sample_id","bin","ocean","depth"]
integrons = integrons.sort_values(by=['sample_id','bin'])[ORF_info + bin_info]
visited = {}

bin_sample_pairs = integrons[bin_info].values.tolist()
ORF_pnps = integrons[ORF_info].values.tolist()

out = []
for i in range(len(bin_sample_pairs)):
	bl, Ol = bin_sample_pairs[i], ORF_pnps[i]
	pari_name = "??".join(bl) # random separator to make sure the has work
	bin_sample_pnps = None
	if pari_name in visited:
		bin_sample_pnps = visited[pari_name]
	else:
		print(f"adding the bin sample pair: {str(bl)}")
		bin_sample_pnps = bin_sample_ranker(sample_id=bl[0], bin=bl[1], ocean=bl[2], depth=bl[3])
		visited[pari_name] = bin_sample_pnps
	rank, total = bin_sample_pnps.find_pnps_ranking(sample_id=bl[0], bin=bl[1], this_pnps=Ol[1])
	out.append(Ol + bl+ [rank, total])

out_df = pd.DataFrame(out)
out_df.columns = ORF_info + bin_info + ["total_pnps", "rank"]
out_df.to_csv(path_or_buf="MAG_integron_pnps_rank.csv", index=False)
