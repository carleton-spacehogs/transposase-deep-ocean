import os


csg = "Major curlin subunit"
lap = "Lipopolysaccharide assembly protein"
pel = "Pellicle polysaccharide"
# lapABCE, csgABCD, pelABCD, tapA, tasA, esp, and FapC 
name_resolver = {"csgA":f"{csg} A", "csgB":f"{csg} B", "csgC":f"{csg} C", "csgD":f"{csg} D", 
"lapA":f"{lap} A", "lapB":f"{lap} B","lapC":f"{lap} C", "lapE":f"{lap} E",
"pelA":f"{pel} A", "pelB":f"{pel} B", "pelC":f"{pel} C", "pelD":f"{pel} D",
"fasA": "Fimbrial_protein 987P", "fapC": "FapC amyloid-like fimbriae protein"}

def resolve_name(gene_short):
	clean_short, full_name = "null", "null"
	if len(gene_short) == 4:
		clean_short = gene_short[:-1].lower() + gene_short[-1].upper()
		return clean_short, name_resolver[clean_short] # full name
	elif gene_short == "ef0056" or gene_short == "esp":
		return "esp", "Enterococcal surface protein"
	elif gene_short == "TAPA_BACSU":
		return "TapA", "TasA anchoring/assembly protein"
	elif gene_short == "TASA_BACSU": 
		return "TasA", "Major biofilm matrix component"
	else:
		return clean_short, full_name

all_genes = []
files = "Biofilm_seeds_forPublication.txt"
with open(files) as genes: 
	all_genes = genes.readlines()

return_list = []
for line in all_genes:
	if line.startswith(">sp|") or line.startswith(">tr|"):
		accession_id, more_info = line.split("|")[1], line.split("|")[2]
		gene_short = (more_info.split("GN=")[-1]).split(" ")[0]
		clean_short, full_name = resolve_name(gene_short)
		return_list.append(["Uniprot", clean_short, full_name, accession_id])
	elif line.startswith(">K"):
		info = line.strip(">").strip("\n").split("_") #>K21007_PelB
		accession_id, gene_short = info[0], info[1]
		return_list.append(["KEGG", clean_short, full_name, accession_id])
	elif line.startswith(">"):
		print("do this manually" + line) 

sorted_output = sorted(return_list, key=lambda r: r[2])
sorted_output.insert(0, ["source", "abbr.", "full gene name", "accession number"])

import csv
with open('biofilm_seeds_accession.csv','w') as f:
	writer = csv.writer(f)
	writer.writerows(sorted_output)

files = "Transposase_Database.faa"
all_genes = []
return_trans = set()
with open(files) as genes: 
	all_genes = genes.readlines()

count = 0
for line in all_genes:
	if line.startswith(">"):
		count += 1
		info = line.strip(">").strip("\n").split("/") # >B9TL16_RICCO/11-65
		return_trans.add(info[0])
print("total entries: " + str(count))

count = 0
accession_per_row = 19 # it is actually 20
with open('pfam_transposase_seeds_accession.csv','w') as f:
	for gene in return_trans:
		count += 1
		f.write(gene + ",")
		if count % accession_per_row == 0: 
			f.write("\n")

