#!/bin/bash
#
#Script to mark duplicates in BAM files using picard
#
for bam_file in ../data/processed_data/bru-mor_bam/*.sorted.bam; do
    bam_file_name=$(basename "$bam_file" .sorted.bam)

		/opt/apps/picard-tools/1.119/bin/MarkDuplicates \
			I=$bam_file \
			O=${bam_file_name}.dupMarked.bam \
			M=${bam_file_name}.dupMetrics.txt
done