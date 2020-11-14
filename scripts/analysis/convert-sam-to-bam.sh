# #!/bin/bash
#
#Script to convert sam files to sorted bam files
#
# USAGE:
# ./convert-sam-to-bam.sh \
# ../../data/scratch/reads-alignment-output/savage \
# ../../data/scratch/sam-to-bam-output/savage \
# ../../data/scratch/reads-alignment-output/telleria/SRR965341.sam \
# ../../data/scratch/sam-to-bam-output/telleria

# make directory for sorted bam files output
mkdir -p ../../data/scratch/sam-to-bam-output/savage
mkdir -p ../../data/scratch/sam-to-bam-output/telleria

# sam files directory
SAM_DIR=$1

# bam files directory
BAM_DIR=$2

# telleria files and output directory
SAM_FILE=$3
OUT_DIR=$4

# convert sam file to sorted bam files
for sam_file in ${SAM_DIR}/*.sam; do
	sam_file_name=$(basename "$sam_file" .sam)
        samtools view -S -b $sam_file | \
        samtools sort -o ${BAM_DIR}/${sam_file_name}.sorted.bam
done

# telleria paired-end data processing
# file name
tel_sam_file_name=$(basename "$SAM_FILE" .sam)

samtools view -S -b $SAM_FILE | \
samtools sort -o ${OUT_DIR}/${tel_sam_file_name}.sorted.bam






