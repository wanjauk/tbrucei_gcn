# #!/bin/bash
#
#Script to generate MultiQC reports using the MultiQC tool
#
# USAGE:
# bash generate-multiqc-reports.sh \
# ../../results/figures/fastqc_reports/savage \
# ../../data/scratch/reads-alignment-output/savage \
# ../../data/intermediate/tbrucei_read_counts/savage \
# ../../results/figures/multiqc_reports/savage

# bash generate-multiqc-reports.sh \
# ../../results/figures/fastqc_reports/telleria \
# ../../data/scratch/reads-alignment-output/telleria \
# ../../data/intermediate/tbrucei_read_counts/telleria \
# ../../results/figures/multiqc_reports/telleria

# make directory to store the results
mkdir -p ../../results/figures/multiqc_reports/savage/ 
mkdir -p ../../results/figures/multiqc_reports/telleria/

# fastqc directory
FASTQC_DIR=$1

# Hisat2 directory
HISAT2_DIR=$2

# Htseq directory
HTSEQ_DIR=$3

# Multiqc output directory
OUT_DIR=$4

# run mutiqc
multiqc ${FASTQC_DIR}/*_fastqc.zip ${HISAT2_DIR}/*.txt ${HTSEQ_DIR}/*.counts.txt --outdir ${OUT_DIR}




