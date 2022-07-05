import csv
import glob

# working dir: /workspace/data/zhongj/Transposase_Project/particle_lifestyle/transposase-deep-ocean/particle-association-redux

# outfile
secretory_CAZ = "./secretory_CAZenzyme.txt"
deep_secretory_CAZ = "./deep_secretory_CAZenzyme.txt"
secretory_pep = "./secretory_peptidase.txt"

CAZ_signalp_outfiles = glob.glob("../../CAZenzyme_signalp_out/*_summary.signalp5")
deep_CAZ_signalp_f = glob.glob("../../deep_metagenome/*_summary.signalp5")
''' a summary.signalp5 output file structure:
# SignalP-5.0   Organism: gram- Timestamp: 20220703215412
# ID    Prediction      SP(Sec/SPI)     TAT(Tat/SPI)    LIPO(Sec/SPII)  OTHER   CS Position
OM-RGC.v2.010671536     LIPO(Sec/SPII)  0.175084        0.000557        0.722455        0.101903        CS pos: 24-25. FCG-LS. Pr: 0.4868
OM-RGC.v2.010591442     OTHER   0.006607        0.001668        0.000785        0.990940
OM-RGC.v2.012952042     OTHER   0.003266        0.000160        0.000483        0.996091
'''

def get_all_signalp(signalp_file):
	signalp = list(csv.reader(open(signalp_file), delimiter='\t'))[2:]
	out = []
	for line in signalp:
		Gene_ID, signal_type = line[0], line[1]
		if signal_type in ['SP(Sec/SPI)','TAT(Tat/SPI)','LIPO(Sec/SPII)']:
			out.append(Gene_ID)
	return out

def main():
	# NOTE!!! this file might contain repeated singal protein ID
	# gramNegative and gramPositive prediction sometimes overlap
	with open(secretory_CAZ, 'w') as a:
		for cf in CAZ_signalp_outfiles:
			CAZ_list = get_all_signalp(cf)
			for Gene_ID in CAZ_list: a.write("%s\n" % Gene_ID)
	
	with open(deep_secretory_CAZ, 'w') as a:
		for cf in deep_CAZ_signalp_f:
			CAZ_list = get_all_signalp(cf)
			for Gene_ID in CAZ_list: a.write("%s\n" % Gene_ID)
		# the same thing will repeat for secretory_peptidase.txt

if __name__ == "__main__":
	main()

