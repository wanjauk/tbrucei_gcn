#!/bin/bash
#
#Script to convert sam files to bam files
#
for sam_file in ~/tbrucei_gcn/data/processed_data/*.sam; do
	sam_file_name=$(echo $sam_file | cut -f1 -d '.')
		samtools view -S -b $sam_file > ${sam_file_name}.bam
done
