import csv
import sys
import pandas as pd

bin_coverages = list(csv.reader(open('../Figure Generating/data/TOBG-Combined.RPKM.csv', 'r')))
master_bin = list(csv.reader(open('../Figure Generating/data/bin_taxon.csv', 'r')))[1:]
bin_of_interested = [row[0] for row in master_bin]

# RPKM (reads per kilobase pair MAG per Million pair metagenome)
read_count_file = './Table1_ReadsPrimaryContigs_modifiedfrom_Tully_et_al_2018.csv'
read_counts = list(csv.reader(open(read_count_file, 'r')))[1:]
read_counts_dict = {row[3]:int(row[4]) for row in read_counts}
cal_depth_reads = pd.read_csv(read_count_file)
depth_sum = cal_depth_reads.groupby(["depth"])["num_reads"].agg('sum')
normaling_factor = [depth_sum["SRF"], depth_sum["DCM"], depth_sum["MES"]]

size_sum = cal_depth_reads.groupby(["size_fraction"])["num_reads"].agg('sum')
# viral (<0.22 μm), girus (0.22–0.8 μm), bacterial (0.22–1.6 μm), and protistan (0.8–5.0 μm)
upper_cut, lower_cut = 2, 0.5
nor_fac_size = [size_sum["prot"], size_sum["bact"]]

samples = bin_coverages[0]
out_file = [["bin", "depth", "size_fraction", "sum_SRF", "sum_DCM", "sum_MES", "sum_particle", "sum_planktonic"]]

par_count, plank_count, SRF_count, DCM_count, MES_count = 0,0,0,0,0

def get_depth(depth_r):
	# depth_r = [RPKM SRF, DCM, MES]
	global SRF_count 
	global DCM_count 
	global MES_count 
	if depth_r[0] + depth_r[1] < depth_r[2]: 
		MES_count += 1
		return "MES"
	elif depth_r[1] + depth_r[2] < depth_r[0]: 
		DCM_count += 1
		return "SRF"
	elif depth_r[0] + depth_r[2] < depth_r[1]: 
		SRF_count += 1
		return "DCM"
	else: 
		return "unsure"

def get_depth2(sum_row_depth, normaling_factor):
	global SRF_count
	global DCM_count
	global MES_count
	def comp(r_sum, nor, ind): # r_sum: read sum
		this_depth_rpkm = r_sum[ind]/nor[ind]
		other_depths_rpkm = sum(r_sum[:ind] + r_sum[ind+1:])/sum(nor[:ind] + nor[ind+1:])
		if this_depth_rpkm/other_depths_rpkm > 2:
			return True
	if comp(sum_row_depth, normaling_factor, 0):
		SRF_count += 1
		return "SRF"
	elif comp(sum_row_depth, normaling_factor, 1): 
		DCM_count += 1
		return "DCM"
	elif comp(sum_row_depth, normaling_factor, 2): 
		MES_count += 1
		return "MES"
	else: 
		return "unsure"


def get_size(size_r):
	global par_count 
	global plank_count 
	# size_r = [particle, planktonic]
	if sum_nor_size[0] / sum_nor_size[1] > upper_cut: 
		par_count += 1
		return "particle"
	elif sum_nor_size[0] / sum_nor_size[1] < lower_cut: 
		plank_count += 1
		return "planktonic"
	else: 
		return "mixed"

# get MAG's depth
for MAG_row in bin_coverages[1:]:
	bin_name = MAG_row[0]
	if bin_name not in bin_of_interested: continue

	out_row = [bin_name]
	sum_row_depth = [0, 0, 0] # SRF, DCM, MES
	sum_row_size = [0, 0] # particle, planktonic
	# print(MAG_row)
	for col in range(1, len(samples)):
		sample = samples[col] # ie. sample = tara151_prot_DCM
		depth, size_type = sample.split("_")[2], sample.split("_")[1]

		coverage = float(MAG_row[col])*read_counts_dict[sample]
		# print(f"cannot find number of reads for {sample}, using average instead")

		if depth == "SRF": sum_row_depth[0] += coverage
		elif depth == "DCM": sum_row_depth[1] += coverage
		elif depth == "MES": sum_row_depth[2] += coverage
		elif depth in ["mixed", "epi"]: pass
		else: sys.exit(f"error! got depth of {depth}")

		if size_type in ["prot"]: sum_row_size[0] += coverage
		elif size_type in ["bact"]: sum_row_size[1] += coverage

	sum_normalized_depth = [a / b for a, b in zip(sum_row_depth, normaling_factor)]
	sum_nor_size = [a / b for a, b in zip(sum_row_size, nor_fac_size)]

	# out_row.append(get_depth(sum_normalized_depth))
	out_row.append(get_depth2(sum_row_depth, normaling_factor))
	out_row.append(get_size(sum_nor_size))

	out_row.extend(sum_normalized_depth)
	out_row.extend(sum_nor_size)
	out_file.append(out_row)

with open("Tara_bins_origin.csv", "w", newline="") as f:
	writer = csv.writer(f)
	writer.writerows(out_file)

print("total:", len(bin_of_interested))
print("SRF:", SRF_count, ", DCM:", DCM_count, ", MES", MES_count)

print("with upper cutoff of:", upper_cut, ", lower cutoff of:", lower_cut)
print(par_count, "MAGs have particle lifestyle,", plank_count, "are plankontic")
