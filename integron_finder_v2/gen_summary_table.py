import sys
import pandas as pd
import numpy as np
from merge_integron_output import deep_rt, anvio_db_dir, oceans

pattern = "nan|Function unknown|General function prediction only"
# merged_integrons[~merged_integrons.function.str.contains(pattern)] # 7966
# merged_integrons[merged_integrons.function.str.contains("Defense mechanisms")] # 734

COG_integrons = pd.read_csv("metagenome-integron-with-COG-functions.csv")
functional_integrons = COG_integrons[~COG_integrons.function.str.contains(pattern)]
total_integron_count = len(functional_integrons.index)

integron_sum = pd.DataFrame(functional_integrons.function.value_counts())
integron_sum['integron_prop'] = integron_sum.function/total_integron_count

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
	whole_func = whole_func[~whole_func.function.str.contains("!!!|" + pattern)]
	all_func = all_func.append(whole_func)

total_all_count = len(all_func.index)
all_sum = pd.DataFrame(all_func.function.value_counts())
all_sum['all_prop'] = all_sum.function/total_all_count

all_sum.columns = ["normal_count", "normal_prop"]
integron_sum.columns = ["integron_count", "integron_prop"]

merged_summary = integron_sum.join(all_sum)
merged_summary = merged_summary.rename_axis('COG_function').reset_index(level=0)
merged_summary.to_csv("integron_known_COG_category_summary.csv", index=False)
