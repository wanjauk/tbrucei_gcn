#!/bin/bash
#
#Script to align reads on the reference genome
#Created April 1, 2019 by Kennedy Mwangi
#
module load star/2.7.0e

#Create output directory if it doesn't exist
mkdir -p ../results/STAR_align_output

GENOME_INDEX=$1 #path to STAR genome index directory.

for fastq in ../data/raw_data/*.fastq; do
    basename=$(echo $fastq | cut -f1 -d '.')

    /opt/apps/star/2.7.0e/bin/STAR \
        --runThreadN 10 \
        --genomeDir $GENOME_INDEX \
        --readFilesIn $fastq \
        --outFileNamePrefix ../results/STAR_align_ouput/${basename}
done
