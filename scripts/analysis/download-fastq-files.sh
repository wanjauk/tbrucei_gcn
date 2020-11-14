#!/bin/bash
#
#Script to download fastq files from European Nucleotide Archive
#
#USAGE: 
# ./download-fastq-files.sh \
# ../../data/raw/savage/savage.fastq.urls.txt /data/kwanjau/savage ;
# ./download-fastq-files.sh \
# ../../data/raw/telleria/telleria.fastq.urls.txt /data/kwanjau/telleria
#
FILE=$1 #File containing fastq url links to EBI FTP site

OUT_DIR=$2 

cat ${FILE} | xargs -n1 wget $3 -P ${OUT_DIR}


#decompress fastq.gz files
#
FASTQ_FILES=${OUT_DIR}/*.fastq.gz

for file in ${FASTQ_FILES}; do
    gunzip ${file}
done
