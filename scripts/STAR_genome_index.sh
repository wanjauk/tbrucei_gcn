#!/bin/bash
#
#Script to create T. brucei genome index using STAR.
#
#Created on March 29, 2019 by Kennedy Mwangi
#
mkdir -p ~/tbrucei_gcn/data/STAR_genome

module load star/2.7.0e

INDEX_DIR=~/tbrucei_gcn/data/STAR_genome/
GENOME_DIR=$1
ANNOTATION_FILE=$2

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir ${INDEX_DIR} \
--genomeFastaFiles ${GENOME_DIR} \
--sjdbGTFfile ${ANNOTATION_FILE} \
--sjdbOverhang 74 \
--sjdbGTFtagExonParentTranscript Parent #GFF file being used instead of GTF
