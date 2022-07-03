import csv
import glob

all_CAZ = "../../diamond_only_dbcan_overview.txt"
secretory_CAZ = "./secretory_CAZenzyme.txt"
secretory_pep = "./secretory_peptidase"

all_CAZ = list(csv.reader(open(all_CAZ), delimiter='\t'))
