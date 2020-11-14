# #!/bin/bash
#
#Script to mark duplicates in BAM files using picard
#
# USAGE:
# ./mark-duplicates.sh \
# ../../data/scratch/sam-to-bam-output/savage \
# ../../data/scratch/mark-duplicates-output/savage

# ./mark-duplicates.sh \
# ../../data/scratch/sam-to-bam-output/telleria \
# ../../data/scratch/mark-duplicates-output/telleria

# crate a directory for mark duplicates output
mkdir -p ../../data/scratch/mark-duplicates-output/savage
mkdir -p ../../data/scratch/mark-duplicates-output/telleria

# sorted bam files directory
SORTED_BAM_DIR=$1

# mark duplicates output
MARK_DUPES_OUT=$2

for bam_file in ${SORTED_BAM_DIR}/*.sorted.bam; do
    bam_file_name=$(basename "$bam_file" .sorted.bam)

        picard MarkDuplicates \
			I=$bam_file \
			O=${MARK_DUPES_OUT}/${bam_file_name}.dupMarked.bam \
			M=${MARK_DUPES_OUT}/${bam_file_name}.dupMetrics.txt
done
