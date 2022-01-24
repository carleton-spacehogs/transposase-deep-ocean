import glob
import os
import sys
import csv
import random

depth=sys.argv[1]
all_blastp=glob.glob(f'Tara*_{depth}_0.22-3/*_anvi-genecall_transposase.blastp')

out=[["gene_caller_id", "pnps", "depth", "gene_type"]]

def target_id_set(blastp_file):
	target_id = set()
	with open(blastp_file, newline='') as f:
		blast_result=list(csv.reader(f,delimiter='\t'))
	for line in blast_result:
		target_id.add(line[1])
	return target_id

def defense_set(cog_file):
	defense_id = set()
	with open(cog_file, newline='') as f:
		all_anvi_cog=list(csv.reader(f,delimiter='\t'))
	for line in all_anvi_cog:
		if "Defense" in line[3]: defense_id.add(line[0])
	return defense_id

for blastp in all_blastp:
	biofilm_blastp=blastp.replace("transposase","biofilm")
	all_cog=blastp.rsplit('/', 1)[0] + "/all-COG-categories.txt"
	transposases = target_id_set(blastp)
	biofilms = target_id_set(biofilm_blastp)
	defenses = defense_set(all_cog)
	print(len(defenses))

	pnps_file=blastp.replace("_anvi-genecall_transposase.blastp","_pn-ps/pNpS.txt")
	with open(pnps_file, newline='') as f:
		pnps_result=list(csv.reader(f,delimiter='\t'))[1:]
	for line in pnps_result:
		gene_caller_id, pnps = line[1], line[3]
		if not pnps: continue
		elif float(pnps) > 10: continue
		elif gene_caller_id in transposases:
			out.append([gene_caller_id, pnps, depth, "transposase"])
		elif gene_caller_id in biofilms:
			out.append([gene_caller_id, pnps, depth, "biofilm"])
		elif gene_caller_id in defenses:
			out.append([gene_caller_id, pnps, depth, "defense"])
		else:
			if random.randint(0, 25) == 1:
				out.append([gene_caller_id, pnps, depth, "n"])

t_count = t_sum = d_count = d_sum =0
normal_sum=0
for line in out[1:]:
	if line[3] == "transposase": 
		t_count += 1
		t_sum += float(line[1])
	elif line[3] == "defense": 
		d_count += 1
		d_sum += float(line[1])
	else:
		normal_sum += float(line[1])

print("there are:", t_count, "transposases, pnps average:", t_sum/t_count)
print("there are:", d_count, "defense mech, pnps average:", d_sum/d_count)
print("normal pnps average:", normal_sum/(len(out)-t_count-d_count), "total: ", len(out))

with open(f'trans_biofilm_defense_and_normal_pnps_{depth}.csv','w') as f:
	writer = csv.writer(f)
	writer.writerows(out)
