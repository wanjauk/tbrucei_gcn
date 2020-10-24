#!/bin/bash
#
#Script to align reads to the indexed genome using HISAT2
#
for fastq in ../data/raw_data/*.fastq; do
    fqname=$(basename "$fastq" .fastq)

		hisat2 \
		 -x bru-mor_genome_index_hisat2 \
		 -U ${fastq} \
		 -S ${fqname}.sam \
		 -p 8 \
		--summary-file ${fqname}.txt \
		--new-summary
done

#move the output sam files to a new directory

#create directory if not exists
mkdir -p ../data/processed_data/bru-mor_sam

#move the sam files
mv ../data/raw_data/*.sam ../data/processed_data/bru-mor_sam/
