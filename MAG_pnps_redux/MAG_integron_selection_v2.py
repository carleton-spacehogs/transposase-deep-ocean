import pandas as pd
import csv
from MAG_singleCopyGene_defense_pnps_ratio import ocean_depths, read_pnps

'''MAG_integron_selection.py was wrong (didn't get the right result)
However, it has helper function that led to the right thing in
../integron_finder_v2/merge_integron_output.py; and
../integron_finder_v2/MAG-integron-all.csv is correct
'''

integrons = list(csv.reader(open("../integron_finder_v2/MAG-integron-all.csv")))
integron_anvi_id = {}
for l in integrons:
	anvi_id, ocean, COG_accession, COG_function = l[0], l[-3], l[-2], l[-1]
	MAG_name = l[1].split('_', 1)[0]
	integron_anvi_id[f"{ocean}-{anvi_id}"] = [MAG_name,COG_accession,COG_function]

MAG_rt="/researchdrive/zhongj2/MAG_pnps_redux"

pnps_integron = []
for ocean, depths in ocean_depths().items():
	for depth in depths:
		if depth != "deep":
			all_pnps = read_pnps(MAG_rt, ocean, depth).values.tolist()
			for l in all_pnps:
				area_id = f"{ocean}-{l[0]}"
				if area_id in integron_anvi_id:
					new_row = l + [ocean, depth] + integron_anvi_id[area_id]
					pnps_integron.append(new_row)

for ocean in ["IN", "SAT", "NAT", "SP", "NP"]:
	all_pnps = read_pnps(MAG_rt, "deep", ocean).values.tolist()
	for l in all_pnps:
		area_id = f"deep-{l[0]}"
		if area_id in integron_anvi_id:
			l[2] = l[2].split("_")[1] # NP_SRR3965585
			new_row = l + [ocean, "deep"] + integron_anvi_id[area_id]
			pnps_integron.append(new_row)

pnps_integron_df = pd.DataFrame(pnps_integron)
pnps_integron_df.columns = ["gene_callers_id","pnps","sample_id","ocean","depth","bin","COG_accession","COG_category"]

all_MAG_info = pd.read_csv(f"MAG_info_db/MAG_all.csv").astype(str)
pnps_integron_MAG = pnps_integron_df.merge(all_MAG_info, on=["bin","sample_id","depth"])
# pnps_integron_MAG[pnps_integron_MAG.total_count.astype(float) > 10]

pnps_integron_MAG.to_csv('MAG_integron_pnps_all_v2.csv', sep=',', index=False)

# count = 0
# for l in pnps_integron:
# 	if "Defense" in l[-1]:
# 		count += 1

# if __name__ == "__main__":
# 	main()
