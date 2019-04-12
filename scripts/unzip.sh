#/bin/bash
#
#Script to decompress fastq.gz files
#Created March 28,2019 by Kennedy Mwangi
#
FASTQ_FILES=../data/raw_data/*.fastq.gz

for file in ${FASTQ_FILES}; do
    gunzip ${file}
done
