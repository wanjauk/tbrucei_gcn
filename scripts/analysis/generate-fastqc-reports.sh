# #!/bin/bash
#
#Script to generate FastQC reports using the FastQC tool
#
# USAGE:
# bash generate-fastqc-reports.sh ../../results/figures/fastqc_reports/savage ../../data/raw/savage 
# bash generate-fastqc-reports.sh ../../results/figures/fastqc_reports/telleria ../../data/raw/telleria
#
# make directory to store the results
mkdir -p ../../results/figures/fastqc_reports/savage/; 
mkdir -p ../../results/figures/fastqc_reports/telleria/
#
#load fastqc module
# module load fastqc/0.11.4

# fastqc reports directory
REPORT_DIR=$1

# fastq files directory
FASTQ_DIR=$2

for file in $FASTQ_DIR/*.fastq; do
   fastqc ${file} -o ${REPORT_DIR} -f fastq
done
