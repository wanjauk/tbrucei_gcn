# #!/bin/bash
#
#Script to align reads to the indexed genome using HISAT2
#

#USAGE:
# ./reads-alignment.sh \
# ../../data/raw/savage \
# ../../data/scratch/indexed_genome \
# ../../data/scratch/reads-alignment-output/savage \
# ../../data/raw/telleria/SRR965341_1.fastq \
# ../../data/raw/telleria/SRR965341_2.fastq \
# ../../data/scratch/reads-alignment-output/telleria

# create alignment output directory
mkdir -p ../../data/scratch/reads-alignment-output/savage
mkdir -p ../../data/scratch/reads-alignment-output/telleria

# Fastq directory
FASTQ_DIR=$1

# genome index directory
INDEX_DIR=$2

# alignment output directory
ALIGN_OUT=$3

# telleria paired-end reads
READ1=$4
READ2=$5
TELLERIA_OUT=$6

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


# Process paired-end data (Telleria's study)
telleria_fastq_name=$(basename "$READ1" _1.fastq)

hisat2 \
		 -x ${INDEX_DIR}/bru-mor_genome_index_hisat2 \
		 -1 ${READ1} \
         -2 ${READ2} \
		 -S ${TELLERIA_OUT}/${telleria_fastq_name}.sam \
		 -p 6 \
		--summary-file ${TELLERIA_OUT}/${telleria_fastq_name}.txt \
		--new-summary

