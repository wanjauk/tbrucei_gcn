#!/bin/bash
#
#Script to run FastQC reports using the FastQC tool
#Written on March 28, 2019 by Kennedy Mwangi
#
#load fastqc module
module load fastqc/0.11.4

FASTQ_DIR=../data/raw_data/*.fastq

#create output directory if it doesn't exist.
mkdir -p ../results/fastqc_reports

REPORTS_DIR=../results/fastqc_reports/

for file in ${FASTQ_DIR}; do
   fastqc -f fastq -o ${REPORTS_DIR} ${file}
done
