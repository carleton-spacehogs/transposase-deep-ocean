import sys
import pandas as pd
import numpy as np
sys.path.insert(0, '../MAG_pnps_redux')
from MAG_integron_selection import explode_COG, integron_dict_helper, merge_anvi_gene_call
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths

deep_rt = "/researchdrive/zhongj2/deep_ocean_bins/"
anvio_db_dir="/researchdrive/zhongj2/Tara_Oceans_Bins"
pattern = "nan|Function unknown|General function prediction only"
o_depth = ocean_depths()
oceans = list(o_depth.keys()) + ["deep"]

COG_integrons = pd.read_csv("metagenome-integron-with-COG-functions.csv").astype("str") # 25178 rows
COG_integrons = COG_integrons[~COG_integrons.function.str.contains(pattern)]
# COG_integrons[COG_integrons.function == "nan"] # 16177 rows, 16177/25178 = 0.6425054
# COG_integrons['function'] = COG_integrons.function.str.replace(pattern, "Function unknown", regex=True)
# COG_integrons[COG_integrons.function == "Function unknown"] # 17439 rows, 17439/25178 = 0.69262848518
integron_sum = pd.DataFrame(COG_integrons.function.value_counts())
integron_sum['integron_prop'] = integron_sum.function/len(COG_integrons.index)

def get_anvi_COG_short(root, ocean):
	COG_file = f"{root}/{ocean}_bins_v2/all-COG-categories.txt"
	if ocean == "deep":
		COG_file = f"{deep_rt}/all-anvi-cog-category.tsv"
	category = pd.read_csv(COG_file, sep="\t", usecols=["function"])
	return category

all_func = pd.DataFrame()
for ocean in oceans:
	print(ocean)
	whole_func = get_anvi_COG_short(anvio_db_dir, ocean)
	all_func_count = len(whole_func.index)
	whole_func = whole_func[~whole_func.function.str.contains("!!!|" + pattern)]
	# print("annotation rate: " + str(len(whole_func.index)/all_func_count))
	all_func = all_func.append(whole_func)

total_all_count = len(all_func.index)
all_sum = pd.DataFrame(all_func.function.value_counts())
all_sum['all_prop'] = all_sum.function/total_all_count

all_sum.columns = ["normal_count", "normal_prop"]
integron_sum.columns = ["integron_count", "integron_prop"]

merged_summary = integron_sum.join(all_sum)
merged_summary = merged_summary.rename_axis('COG_function').reset_index(level=0)
merged_summary.to_csv("integron_known_COG_category_summary.csv", index=False)
# merged_summary.to_csv("integron_all_COG_category_summary.csv", index=False)
