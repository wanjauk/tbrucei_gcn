#!/bin/bash
#
#Script to index genome using HISAT2
#
GENOME_FILE=$1

hisat2-build ${GENOME_FILE} bru-mor_genome_index_hisat2
