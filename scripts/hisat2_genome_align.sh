#!/bin/bash
#
#Script to align T. brucei reads to the indexed genome using HISAT2
#created April 9, 2019
#
module load hisat/2-2.1.0

for fastq in $(ls ~/tbrucei_gcn/data/raw_data/*.fastq); do
    basename=$(echo $fastq | cut -f1 -d '.')
        hisat2 \
            -x tbrucei_genome_index_hisat2 \
            -U ${fastq} \
            -S ~/tbrucei_gcn/results/hisat2_align_results/${basename}.sam \
            -p 8 \
            --new-summary \

