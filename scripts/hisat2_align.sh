#!/bin/bash
#
#Script to align T. brucei reads to the indexed genome using HISAT2
#
module load hisat/2-2.1.0

#change directory to that of the indexed genome
cd ../data/HISAT2_indexed_genome/

for fastq in ../raw_data/*.fastq; do
    fqname=$(echo $fastq | cut -f1 -d '.')

		hisat2 \
		 -x tbrucei_genome_index_hisat2 \
		 -U ${fastq} \
		 -S ${fqname}.sam \
		 -p 8 \
		--summary-file ${fqname}.txt \
		--new-summary
done

#make directories and move created files into them

mkdir -p ../processed_data
mkdir -p ../../results/hisat2_alignment_summary

mv ../raw_data/*.sam ../processed_data/
mv ../raw_data/SRR*.txt ../../results/hisat2_alignment_summary/
