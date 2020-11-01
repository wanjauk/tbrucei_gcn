# #!/bin/bash
#
#Script to convert sam files to sorted bam files
#
# USAGE:
# ./convert-sam-to-bam.sh \
# ../../data/scratch/reads-alignment-ouput/savage \
# ../../data/scratch/sam-to-bam-output/savage

# make directory for sorted bam files output
mkdir -p ../../data/scratch/sam-to-bam-output/savage

# sam files directory
SAM_DIR=$1

# bam files directory
BAM_DIR=$2

# convert sam file to sorted bam files
for sam_file in ${SAM_DIR}/*.sam; do
	sam_file_name=$(basename "$sam_file" .sam)
# 		samtools view -S -b $sam_file > ${sam_file_name}.bam
        samtools view -S -b $sam_file | \
        samtools sort -o ${BAM_DIR}/${sam_file_name}.sorted.bam
done

# # sort BAM file
# #
# for bam_file in ../data/processed_data/bru-mor_bam/*.bam; do
#     bam_file_name=$(basename "$bam_file" .bam)
# 		samtools sort $bam_file -o ${bam_file_name}.sorted.bam
# done
