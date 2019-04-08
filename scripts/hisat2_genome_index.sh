#!/bin/bash
#
#Script to index T. brucei genome using HISAT2
#Created April 8, 2019 by Kennedy Mwangi
#
module load hisat/2-2.1.0

GENOME_FILE=$1


hisat2-build ${GENOME_FILE} tbrucei_genome_index_hisat2
