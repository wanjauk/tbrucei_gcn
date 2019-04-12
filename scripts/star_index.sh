#!/bin/bash
#
#Script to create T. brucei genome index using STAR.
#
#Created on March 29, 2019 by Kennedy Mwangi
#
#Create required genome directory if not exists.
mkdir -p ../data/STAR_genome

module load star/2.7.0e

INDEX_DIR=../data/STAR_genome/
GENOME_FILE=$1
ANNOTATION_FILE=$2

/opt/apps/star/2.7.0e/bin/STAR \
    --runThreadN 10 \
    --runMode genomeGenerate \
    --genomeDir ${INDEX_DIR} \
    --genomeFastaFiles ${GENOME_FILE} \
    --sjdbGTFfile ${ANNOTATION_FILE} \
    --sjdbOverhang 74 \
    --sjdbGTFtagExonParentTranscript Parent
