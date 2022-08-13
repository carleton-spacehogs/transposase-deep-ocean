import pandas as pd
import csv
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, read_pnps

'''MAG_integron_selection.py was wrong (didn't get the right result)
However, it has helper function that led to the right thing in
../integron_finder_v2/merge_integron_output.py; and
../integron_finder_v2/MAG-integron-all.csv is correct
'''

integrons = list(csv.reader(open("../integron_finder_v2/MAG-integron-all.csv")))
# ARS-28919 l[8]-l[0]
# l[-1] is the function
integron_anvi_id = {f"{l[8]}-{l[0]}":l[-1] for l in integrons}
MAG_rt="/researchdrive/zhongj2/MAG_pnps_redux"

pnps_integron = []
for ocean, depths in ocean_depths().items():
	for depth in depths:
		if depth != "deep":
			all_pnps = read_pnps(MAG_rt, ocean, depth).values.tolist()
			for l in all_pnps:
				area_id = f"{ocean}-{l[0]}"
				if area_id in integron_anvi_id:
					new_row = l + [ocean, depth, integron_anvi_id[area_id]]
					pnps_integron.append(new_row)

for ocean in ["IN", "SAT", "NAT", "SP", "NP"]:
	all_pnps = read_pnps(MAG_rt, "deep", ocean).values.tolist()
	for l in all_pnps:
		area_id = f"{ocean}-{l[0]}"
		if area_id in integron_anvi_id:
			new_row = l + [ocean, "deep", integron_anvi_id[area_id]]
			pnps_integron.append(new_row)

pnps_integron_df = pd.DataFrame(pnps_integron)
pnps_integron_df.columns = ["anvi-id", "pnps", "sample_id", "ocean", "depth", "COG_category"]
pnps_integron_df.to_csv('MAG_integron_pnps_all_v2.csv', sep=',', index=False)

# count = 0
# for l in pnps_integron:
# 	if "Defense" in l[-1]:
# 		count += 1

# if __name__ == "__main__":
# 	main()


