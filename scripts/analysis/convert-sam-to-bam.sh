#!/bin/bash
#
#Script to convert sam files to bam files
#
for sam_file in ../data/processed_data/bru-mor_sam/*.sam; do
	sam_file_name=$(basename "$sam_file" .sam)
		samtools view -S -b $sam_file > ${sam_file_name}.bam
done

# move the created bam files to a new directory
mkdir -p ../data/processed_data/bru-mor_bam

mv ../data/processed_data/bru-mor_sam/*.bam ../data/processed_data/bru-mor_bam/


# sort BAM file
#
for bam_file in ../data/processed_data/bru-mor_bam/*.bam; do
    bam_file_name=$(basename "$bam_file" .bam)
		samtools sort $bam_file -o ${bam_file_name}.sorted.bam
done
