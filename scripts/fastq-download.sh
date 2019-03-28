#!/bin/bash
#
#Script to download fastq files from European Nucleotide Archive
#Created March 28,2019 by Kennedy Mwangi
#
FILE=$1 #link to the file containing fastq url links to ENA FTP database
OUTDIR=~/tbrucei_gcn/data/raw_data/

cat ${FILE} | xargs -n1 wget$2 -P ${OUTDIR}
