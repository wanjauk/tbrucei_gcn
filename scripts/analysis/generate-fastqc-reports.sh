#!/bin/bash
#
#Script to run FastQC reports using the FastQC tool
# ./generate-fastqc-reports.sh /data/kwanjau/savage/*.fastq ../../results/figures/fastqc_reports/savage/; \ 
# ./generate-fastqc-reports.sh /data/kwanjau/telleria/*.fastq ../../results/figures/fastqc_reports/telleria/
#
#load fastqc module
# module load fastqc/0.11.4

# fastq files directory
FASTQ_DIR=$1

# fastqc reports directory
REPORTS_DIR=$2

for file in ${FASTQ_DIR}; do
   fastqc -f fastq -o ${REPORTS_DIR} ${file}
done
