# #!/bin/bash

# USAGE: 
# ./concatenate-and-index-genome.sh \
# ../../data/scratch/concatenated_genomes/brucei-morsitans_genomes.fasta \
# ../../data/scratch/indexed_genome

# make a directory to store the concatenated and indexed genomes
 mkdir -p ../../data/scratch/concatenated_genomes
 mkdir -p ../../data/scratch/indexed_genome

# concatenate the genome files
cat ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927_Genome.fasta \
../../data/scratch/glossina/Glossina-morsitans-Yale_SCAFFOLDS_GmorY1.fa \
> ../../data/scratch/concatenated_genomes/brucei-morsitans_genomes.fasta

# concatenate the gtf annotation files
cat ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.gtf \
../../data/scratch/glossina/Glossina-morsitans-Yale_BASEFEATURES_GmorY1.9.gtf \
> ../../data/scratch/concatenated_genomes/brucei-morsitans_annotations.gtf

#
#index genome using HISAT2
#
GENOME_FILE=$1

INDEX_DIR=$2

hisat2-build ${GENOME_FILE} ${INDEX_DIR}/bru-mor_genome_index_hisat2
