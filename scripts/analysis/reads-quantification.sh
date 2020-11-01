# #!/bin/bash
#
#Script to counts the number of reads aligned to T. brucei genome using HTSeq.
#Resource: HTSeq documentation https://htseq.readthedocs.io/en/latest/count.html
#
# module load htseq/0.11.2

# USAGE:
# ./reads-quantification.sh \
# ../../data/scratch/sam-to-bam-output/savage \
# ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.gff \
# ../../data/intermediate/tbrucei_read_counts/savage


# make directory for reads count output from htseq
mkdir -p ../../data/intermediate/tbrucei_read_counts/savage

# path to bam files
BAM_DIR=$1

# T. brucei genome annotation file.
GFF_FILE=$2

# read counts path
READ_COUNTS_DIR=$3

# reads counting
for bam_file in ${BAM_DIR}/*.sorted.bam; do
    bam_file_name=$(basename "$bam_file" .sorted.bam)
    
        htseq-count \
            -f bam \
            -s no \
            -t exon \
            -i Parent \
            -n 6 \
            $bam_file \
            $GFF_FILE \
            > ${READ_COUNTS_DIR}/${bam_file_name}.counts.txt
done
