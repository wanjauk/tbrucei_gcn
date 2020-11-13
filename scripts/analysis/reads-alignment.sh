# #!/bin/bash
#
#Script to align reads to the indexed genome using HISAT2
#

#USAGE:
# ./reads-alignment.sh \
# ../../data/raw/savage \
# ../../data/scratch/indexed_genome \
# ../../data/scratch/reads-alignment-output/savage

# create alignment output directory
mkdir -p ../../data/scratch/reads-alignment-output/savage

# Fastq directory
FASTQ_DIR=$1

# genome index directory
INDEX_DIR=$2

# alignment output directory
ALIGN_OUT=$3

for fastq in ${FASTQ_DIR}/*.fastq; do
    fqname=$(basename "$fastq" .fastq)

		hisat2 \
		 -x ${INDEX_DIR}/bru-mor_genome_index_hisat2 \
		 -U ${fastq} \
		 -S ${ALIGN_OUT}/${fqname}.sam \
		 -p 6 \
		--summary-file ${ALIGN_OUT}/${fqname}.txt \
		--new-summary
done
