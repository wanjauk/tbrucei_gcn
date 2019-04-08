#!/bin/bash
#
#Script to align T. brucei reads to the indexed genome using HISAT2
#created April 9, 2019
#
module load hisat/2-2.1.0

for fastq in  ~/tbrucei_gcn/data/raw_data/*.fastq; do
    fqname=$(echo $fastq | cut -f1 -d '.')
        hisat2 \
            -x tbrucei_genome_index_hisat2 \
            -U ${fastq} \
            -S ${fqname}.sam \
            -p 8 \
            --new-summary ${fqname}.txt
done
