#!/bin/bash
#
#Script to convert sam files to bam files
#
for sam_file in ../data/processed_data/bru-mor_sam/*.sam; do
	sam_file_name=$(echo $sam_file | cut -f1 -d '.')
		samtools view -S -b $sam_file > ${sam_file_name}.bam
done

# move the created bam files to a new directory
mkdir -p ../data/processed_data/bru-mor_bam

mv ../data/processed_data/bru-mor_sam/*.bam ../data/processed_data/bru-mor_bam/
