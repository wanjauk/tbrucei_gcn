#!/bin/bash

# make a directory to store the concatenated genomes
 mkdir -p ../data/brucei-morsitans

# copy the genome files to the created directory and concatenate them
cp ../data/glossina_genome_scaffolds/glossina-* ../data/brucei-morsitans/
cp ../data/tbrucei_genome/*.fasta ../data/brucei-morsitans/
cat ../data/brucei-morsitans/*.fa* > ../data/brucei-morsitans/brucei-morsitans_genomes.fasta

#
#index genome using HISAT2
#
GENOME_FILE=$1

hisat2-build ${GENOME_FILE} bru-mor_genome_index_hisat2
