#!/bin/bash
#
#Script to counts the number of reads aligned to T. brucei genome using HTSeq.
#Resource: HTSeq documentation https://htseq.readthedocs.io/en/latest/count.html
#Created on April 2, 2019 by Kennedy Mwangi
#
module load htseq/0.11.2

mkdir -p ../results/HTSeq_count_results

GFF_FILE=$1

for samfile in $(ls ../results/STAR_align_output/*.??); do
    basename=$(echo $samfile | cut -f1 -d '.') #check the naming of .sam files
    
    htseq \
        --format=sam \
        --stranded=?? \
        --type=exon \ #feature type
        --idattr=Parent \
        --samout=$basename.count \
        $samfile \
        ${GFF_FILE}
done
