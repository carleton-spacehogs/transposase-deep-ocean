import glob
import os
import sys
import csv
import random

depth=sys.argv[1]
all_blastp=glob.glob(f'f_Tara*_{depth}_0.22-3/*_anvi-genecall_transposase.blastp')

out=[["gene_caller_id", "pnps", "gene_type"]]

for blastp in all_blastp:
	transposases = set()
	with open(blastp, newline='') as f:
		blast_result=list(csv.reader(f,delimiter='\t'))
	for line in blast_result:
		transposases.add(line[1])
	pnps_file=blastp.replace("_anvi-genecall_transposase.blastp","_pn-ps/pNpS.txt")
	with open(pnps_file, newline='') as f:
		pnps_result=list(csv.reader(f,delimiter='\t'))[1:]
	for line in pnps_result:
		gene_caller_id, pnps = line[1], line[3]
		if not pnps: continue
		elif float(pnps) > 10: continue
		elif gene_caller_id in transposases:
			out.append([gene_caller_id, pnps, f"{depth}_trans"])
		else:
			if random.randint(0, 50) == 1:
				out.append([gene_caller_id, pnps, depth])

trans_count=0
trans_sum=0
normal_sum=0
for line in out[1:]:
	if line[2] == f"{depth}_trans": 
		trans_count += 1
		trans_sum += float(line[1])
	else:
		normal_sum += float(line[1])

print("there are:", trans_count, "transposases, pnps average:", trans_sum/trans_count)
print("normal pnps average:", normal_sum/(len(out)-trans_count))

with open(f'trans_and_normal_pnps_{depth}.csv','w') as f:
	writer = csv.writer(f)
	writer.writerows(out)
