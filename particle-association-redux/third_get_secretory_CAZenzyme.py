import csv
import glob
import os

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

# outfile
secretory_CAZ = "./secretory_CAZenzyme.txt"
deep_secretory_CAZ = "./deep_secretory_CAZenzyme.txt"
secretory_pep = "./secretory_peptidase.txt"
deep_secretory_pep = "./deep_secretory_peptidase.txt"

CAZ_signalp_outfiles = glob.glob("../../OM-RGC_CAZenzyme_signalp/*_summary.signalp5")
pep_signalp_outfiles = glob.glob("../../OM-RGC_peptidase_signalp/*_summary.signalp5")

deep_CAZ_signalp_f = glob.glob("../../deep_metagenome/deep_CAZenzyme_*_summary.signalp5")
deep_pep_signalp_f = glob.glob("../../deep_metagenome/deep_peptidase_*_summary.signalp5")
''' a summary.signalp5 output file structure:
# SignalP-5.0   Organism: gram- Timestamp: 20220703215412
# ID    Prediction      SP(Sec/SPI)     TAT(Tat/SPI)    LIPO(Sec/SPII)  OTHER   CS Position
OM-RGC.v2.010671536     LIPO(Sec/SPII)  0.175084        0.000557        0.722455        0.101903        CS pos: 24-25. FCG-LS. Pr: 0.4868
OM-RGC.v2.010591442     OTHER   0.006607        0.001668        0.000785        0.990940
OM-RGC.v2.012952042     OTHER   0.003266        0.000160        0.000483        0.996091
'''

def get_all_signalp(signalp_file):
	out = []
	if os.path.getsize(signalp_file) == 0: # the files I created manually,
		# this MAG doesn't have any blastp match to peptidase
		return out
	signalp = list(csv.reader(open(signalp_file), delimiter='\t'))[2:]
	for line in signalp:
		Gene_ID, signal_type = line[0], line[1]
		if signal_type in ['SP(Sec/SPI)','TAT(Tat/SPI)','LIPO(Sec/SPII)']:
			out.append(Gene_ID)
	return out

def get_secretory_file(outfile, signalp_files):
	out = []
	for piece in signalp_files:
		out = out + get_all_signalp(piece)
	# gramNegative and gramPositive prediction sometimes overlap
	out = set(out)
	with open(outfile, 'w') as a:
		[a.write("%s\n" % Gene_ID) for Gene_ID in out]
	print("output: " + outfile)

def main():
	get_secretory_file(secretory_CAZ, CAZ_signalp_outfiles)
	get_secretory_file(secretory_pep, pep_signalp_outfiles)
	
	# deep ocean
	get_secretory_file(deep_secretory_CAZ, deep_CAZ_signalp_f)
	get_secretory_file(deep_secretory_pep, deep_pep_signalp_f)

if __name__ == "__main__":
	main()

