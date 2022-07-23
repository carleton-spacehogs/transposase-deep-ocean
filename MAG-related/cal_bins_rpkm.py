import os
import glob
import csv
# RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
# = total read base count in a MAG / ( geneLength/1000 * total read bases/1,000,000 )

bin_size = []
with open("bin_size_bp.txt") as file:
	bin_size = list(csv.reader(file, delimiter="\t"))[1:]

for r in range(len(bin_size)):
	bin_size[r].append(0.0) # particle abundance
	bin_size[r].append(0.0) # free-living abundance
	bin_size[r].append("NA") # abundance ratio: particle/free-living
	bin_size[r].append("NA") # lifestyle: [if abundance ratio > 2: particle-associate]
										# [if abundance ratio < 0.5: free-living]

total_reads_dict = {}
with open("total_sample_bp.txt") as file:
	lines = list(csv.reader(file, delimiter="\t"))[1:]
	for l in lines:
		total_reads_dict[l[0]] = float(l[1]) * 1000000000

cov_files = glob.glob('../*_bin_cov.txt') # ../SRR3965592_0.2_bin_cov.txt
header = ["bin", "bin_size", "particle_abundance", "freeLiving_abundance", "ratio", "lifestyle"]
for cov_f in cov_files:
	bin_contig_dict = {}
	col_name = cov_f.strip("../").strip("_bin_cov.txt")
	header.append(col_name)
	sample_name = col_name.split("_")[0]
	with open(cov_f, "r") as f:
		contig_covs = f.readlines()
	
	particle = True if "0.8" in col_name else False

	MAG_cov = {}
	for line in contig_covs:
		s = line.split("\t")
		k_name, contig_len, MAG, cov = s[0], s[2], s[3], int(s[4]) # total number base pairs mapped to this contig
		# samtools bedcov: Reports the total read base count (i.e. the sum of per base read depths) for each genomic region specified in the supplied BED file. 
		if MAG in MAG_cov:
			MAG_cov[MAG] += cov
		else:
			MAG_cov[MAG] = cov
	
	total_read_this_sample = total_reads_dict[sample_name]
	for r in range(len(bin_size)):
		MAG, MAG_size = bin_size[r][0], float(bin_size[r][1])
		MAG_cov_this_sample = MAG_cov[MAG]
		RPKM = MAG_cov_this_sample/((MAG_size/1000) * (total_read_this_sample/1000000)) # total read base count in a MAG / ( geneLength/1000 * total read bases/1,000,000 )
		bin_size[r].append(str(RPKM))
		if particle:
			bin_size[r][2] += RPKM
		else:
			bin_size[r][3] += RPKM

final_out = [header]
for row in bin_size:
	abundance_ratio = row[2]/row[3]
	row[4] = str(abundance_ratio) # particle/free-living
	if (abundance_ratio > 2): 
		row[5] = "particle"
	elif (abundance_ratio < 0.5):
		row[5] = "planktonic"
	else:
		row[5] = "mixed"
	row[2] = str(row[2]) # free-living abundance
	row[3] = str(row[3]) # abundance ratio
	final_out.append(row)

outfile = open('Malaspina-Combined.RPKM.csv', 'w')
for line in final_out:
	outfile.write(",".join(line)+"\n")
outfile.close()



