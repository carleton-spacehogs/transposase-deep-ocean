#!/bin/bash
all_db=$(ls *.db)

for db in $all_db; do
	site=${db%".db"}
	if [ ! -d "$site" ]; then
		if [ -f "$site.bam.bai" ]; then
			echo $site
			bamprofile=${site}.bam-ANVIO_PROFILE
			bamdb=${bamprofile}/PROFILE.db
			scv=${site}_SCV.txt
			contigdb=${site}.db
			/Accounts/zhongj2/github/anvio/bin/anvi-profile -i ${site}.bam -c $contigdb --profile-SCVs -T 10 --output-dir $bamprofile
			anvi-script-add-default-collection -p $bamdb -c $contigdb -C DEFAULT -b EVERYTHING
			# /Accounts/zhongj2/github/anvio/bin/anvi-migrate $bamdb --migrate-dbs-safely
			/Accounts/zhongj2/github/anvio/bin/anvi-gen-variability-profile -p $bamdb -c $contigdb --engine CDN -o $scv -C DEFAULT -b EVERYTHING
			anvi-get-pn-ps-ratio -V $scv -c $contigdb -o ${site}_pn-ps -m 20
			mkdir f_${site}
			mv ${site}* ./f_${site}/.
		fi
	fi
done



