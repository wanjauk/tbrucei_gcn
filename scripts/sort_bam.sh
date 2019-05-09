#!/bin/bash
#
#Script to sort BAM file
#
for bam_file in ~/tbrucei_gcn/data/processed_data/bru-mor_bam/*.bam; do
    bam_file_name=$(echo $bam_file | cut -f1 -d '.')
		samtools sort $bam_file -o ${bam_file_name}.sorted.bam
done
