# python3
# exec(open('pnps_integron_transposase_merge_master.py').read())

with open("all_ocean_merged_integrons.csv") as file: hits = file.readlines()[1:]

integrons_set = set((line.split(",")[-1].strip()+ line.split(",")[0]) for line in hits) #ie, ARS125434

with open("../Tara_Oceans_Bins/all_transposases_calls_anvio.csv") as file: trans = file.readlines()[1:]

trans_set = set((line.split(",")[0]+ line.split(",")[1]) for line in trans) #ie, ARS125434

def for_others():
	with open("all_pNpS.csv") as file: all_pNpS = file.readlines()[1:]
	out=[]
	for line in all_pNpS:
		key=line.split(",")[1]+ line.split(",")[2].strip()
		iden1, iden2 = ",N", ",N"
		if key in integrons_set: iden1 = ",Y"
		if key in trans_set: iden2 = ",Y"
		if "inf" not in line: out.append(line.strip() + iden1 + iden2)
	with open("pNpS_except_SPSAT.csv", 'w') as f:f.write("\n".join(out))
	return "pNpS_except_SPSAT.csv"

def for_SP_and_SAT(ocean):
	def formatter(is_integron): 
		'''not all pnps in the integron-pn-ps are belonged to integrons
		many just happens to be on the same contig as the contig with integrons'''
		filename = f"../Tara_Oceans_Bins/{ocean}_bins_v2/integron-pn-ps/pNpS.txt" if is_integron else f"../Tara_Oceans_Bins/{ocean}_bins_v2/split_cal_pNpS/random_300000-pn-ps/pNpS.txt"
		with open(filename) as file: data = file.readlines()[1:]
		extra_out=[]
		for line in data: 
			splitted = line.split("\t")
			gene_call, pnps= splitted[1], splitted[3].strip()
			is_integron_str = "Y" if (ocean + gene_call) in integrons_set else "N"
			is_trans = "Y" if (ocean + gene_call) in trans_set else "N"
			if "inf" not in pnps and len(pnps)>1: 
				extra_out.append(f"{pnps},{ocean},{gene_call},{is_integron_str},{is_trans}") # check it's not a integron already, 0.25857947857086383,ARS,32,N,N
		return extra_out
	integron=formatter(True)
	non_integron=formatter(False)
	with open(f"pNpS_{ocean}.csv", 'w') as f:f.write("\n".join(integron+non_integron))
	return f"pNpS_{ocean}.csv"

import pandas as pd
import numpy as np
import scipy.stats

def load_pandas(filename, sample_size):
	df = pd.read_csv(filename, names = ["pnps", "ocean", "gene-callers-id", "integron", "transposase"])
	# Q3 = np.quantile(df['pnps'], 0.75)
	# Q1 = np.quantile(df['pnps'], 0.25)
	# IQR = Q3 - Q1
	# lower_range = Q1 - 1.5 * IQR
	# upper_range = Q3 + 1.5 * IQR
	# print(f"{lower_range}, {upper_range}, {lower_range + upper_range}")
	# df = df.query("pnps < @upper_range & pnps > @lower_range")
	df = df.query('pnps < 10')
	integron = df[df['integron']=='Y']
	trans = df[df['transposase']=='Y']
	non_trans_integron = df.query('transposase == "N" & integron == "N"').sample(n=sample_size, random_state=119)
	return integron, trans, non_trans_integron



# integron1, trans1, non_trans_integron1 = load_pandas(for_others(), 220000)
# integron2, trans2, non_trans_integron2 = load_pandas(for_SP_and_SAT("SP"), 100000)
# integron3, trans3, non_trans_integron3 = load_pandas(for_SP_and_SAT("SAT"), 80000)
integron1, trans1, non_trans_integron1 = load_pandas("pNpS_except_SPSAT.csv", 220000)
integron2, trans2, non_trans_integron2 = load_pandas("pNpS_SP.csv", 100000)
integron3, trans3, non_trans_integron3 = load_pandas("pNpS_SAT.csv", 80000)

integron = integron1.append(integron2).append(integron3)
trans = trans1.append(trans2).append(trans3)
non_trans_integron = non_trans_integron1.append(non_trans_integron2).append(non_trans_integron3)

print("integron vs non_trans_integron:")
print(scipy.stats.ttest_ind(integron['pnps'], non_trans_integron['pnps']))
print("transposase vs non_trans_integron:")
print(scipy.stats.ttest_ind(trans['pnps'], non_trans_integron['pnps']))

print(integron.pnps.describe())
print(trans.pnps.describe())
print(non_trans_integron.pnps.describe())

total = integron.append(trans).append(non_trans_integron)
total.to_csv('pNpS_total.csv', sep=',', index=False)

# for Tara Overall pnps:
'''
>>> total.pnps.describe()
count    402416.000000
mean          0.192987
std           0.406200
min           0.000000
25%           0.055763
50%           0.107060
75%           0.195060
max           9.949393
Name: pnps, dtype: float64
'''

# for Malaspina pnps:
df = pd.read_csv("pNpS.txt", sep='\t')
df = df.query('pNpS_gene_reference < 5')
df.pNpS_gene_reference.sample(n=400000, random_state=119).describe()
'''
>>> df.pNpS_gene_reference.sample(n=400000, random_state=119).describe()
count    400000.000000
mean          0.811466
std           1.041871
min           0.000000
25%           0.213097
50%           0.455489
75%           1.012134
max           9.999615
Name: pNpS_gene_reference, dtype: float64
'''



