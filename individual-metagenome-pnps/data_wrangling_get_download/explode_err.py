import csv
import sys

rows = list(csv.reader(open('./all_errs.csv', 'r')))[1:372]
out = []
for row in rows:
	ERRs = row[1].split("|")
	if len(row) < 4: continue
	tara_station, depth, size = row[3], row[4], row[5]
	while len(tara_station) < 3:
		tara_station = "0" + tara_station
	meta_info = f"Tara{tara_station}_{depth}_{size}" # 34 SRF 0.1-0.22
	for err in ERRs:
		out.append([err, depth, size, meta_info])

with open("err_exploded.csv", "w", newline="") as f:
	writer = csv.writer(f)
	writer.writerows(out)