#!/bin/bash
#
#Script to align reads to the indexed genome using HISAT2
#
for fastq in ~/tbrucei_gcn/data/raw_data/*.fastq; do
    fqname=$(echo $fastq | cut -f1 -d '.')

		hisat2 \
		 -x bru-mor_genome_index_hisat2 \
		 -U ${fastq} \
		 -S ${fqname}.sam \
		 -p 8 \
		--summary-file ${fqname}.txt \
		--new-summary
done
