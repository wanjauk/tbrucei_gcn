#!/usr/bin/bash
#
#Script to decompress fastq.gz files
#Created March 28,2019 by Kennedy Mwangi
#
DIR=~/tbrucei_gcn/data/raw_data/*.fastq.gz

for file in ${DIR}; do
    gunzip ${file}
done
