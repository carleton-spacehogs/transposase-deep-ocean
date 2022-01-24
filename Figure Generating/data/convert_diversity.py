import skbio

index_dict = {}
index=0
for char in "CEFGHIP":
	index_dict[char]=index
	index += 1

# P, inorganic ion transport and metabolism; 
# C, energy production and conversion; 
# G, carbohydrate metabolism and transport; 
# E, amino acid metabolism and transport; 
# F, nucleotide metabolism and transport; 
# H, coenzyme metabolism; 
# I, lipid metabolism; 

# or C	E	F	G	H	I	P

out = ["bin,COG_group,total_metabolic_ORFs,simpson,shannon_base2", "energy production and conversion",
"amino acid", "nucleotide", "carbohydrate",
"coenzyme", "lipid", "inorganic ion"]
out = [",".join(out)]
with open("tara_bin_COG_diversity.csv") as file: 
	bins = file.readlines()[1:]

for a_bin in bins:
	meta_distribution = a_bin.split(",")[1]
	raw_counts = meta_distribution.split(" ")
	#print(raw_counts)
	counts = [] 
	out_counts = ["0"]*7 # CEFGHIP
	for i in range(len(raw_counts)):
		if i%2 == 1: 
			meta_count = int(raw_counts[i])
			counts.append(meta_count)
			put_to = index_dict[raw_counts[i+1]]
			out_counts[put_to] = str(meta_count)
	simpson_diversity = skbio.diversity.alpha.simpson(counts)
	shannon_diversity = skbio.diversity.alpha.shannon(counts, base=2)
	line = f"{a_bin.strip()},{str(simpson_diversity)},{str(shannon_diversity)},{','.join(out_counts)}" 
	out.append(line)

with open(f"tara_bin_COG_diversity2.csv", 'w') as f:
	f.write("\n".join(out))



