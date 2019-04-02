#!/bin/bash
#
#Script to align reads on the reference genome
#Created April 1, 2019 by Kennedy Mwangi
#
module load star/2.7.0e

mkdir -p ~/tbrucei_gcn/results/STAR_align_output

GENOME_INDEX=$1 #path to genome index directory.

for fastq in $(ls ~/tbrucei_gcn/data/raw_data/*.fastq); do
    basename=$(echo $fastq | cut -f1 -d '.')

    STAR \
        --runThreadN 4 \
        --genomeDir $GENOME_INDEX \
        --readFilesIn $fastq \
        --outFileNamePrefix ../results/STAR_align_ouput/$basename.
