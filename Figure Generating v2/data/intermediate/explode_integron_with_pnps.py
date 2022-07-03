import glob
import os
import csv

integrons = list(csv.reader(open("all_integron_with_pnps.csv")))
out = []
for i in range(len(integrons)):
	func_category = integrons[i][3]
	if func_category and "!!!" in func_category:
		func = integrons[i][2]
		ocean, gene_callers_id, pnps, from_integron_type, log_pnps, gene_type = integrons[i][0], integrons[i][1], integrons[i][4], integrons[i][5], integrons[i][6], integrons[i][7]
		func_category_list = func_category.split("!!!")
		for i in range(len(func_category_list)):
			individual_func = func_category_list[i]
			new_gene_type = "defense" if individual_func == "Defense mechanisms" else "non_defense"
			newline = [ocean, gene_callers_id, func, individual_func, pnps, from_integron_type, log_pnps, new_gene_type]
			out.append(newline)
	else:
		out.append(integrons[i]) # don't do anything

with open("all_integron_with_pnps_exploded.csv", 'w', newline='') as f:
	write = csv.writer(f)
	write.writerows(out)
