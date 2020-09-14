#!/bin/bash
#
#Script to sort BAM file
#
for bam_file in ../data/processed_data/bru-mor_bam/*.bam; do
    bam_file_name=$(basename "$bam_file" .bam)
		samtools sort $bam_file -o ${bam_file_name}.sorted.bam
done
