#!/bin/bash
#
#Script to decompress fastq.gz files
#
FASTQ_FILES=../data/raw_data/*.fastq.gz

for file in ${FASTQ_FILES}; do
    gunzip ${file}
done
