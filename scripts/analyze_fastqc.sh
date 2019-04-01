#/bin/bash
#
#Script to run FastQC reports using the FastQC too
#Written on March 28, 2019 by Kennedy Mwangi
#
DATA_DIR=~/tbrucei_gcn/data/raw_data/*.fastq
mkdir -p ~/tbrucei_gcn/results/fastqc_reports #create output directory if it doesn't exist.

REPORTS_DIR=~/tbrucei_gcn/results/fastqc_reports/

for file in ${DATA_DIR}; do
   fastqc -f fastq -o ${REPORTS_DIR} ${file}
done
