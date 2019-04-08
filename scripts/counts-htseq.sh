#!/bin/bash
#
#Script to counts the number of reads aligned to T. brucei genome using HTSeq.
#Resource: HTSeq documentation https://htseq.readthedocs.io/en/latest/count.html
#Created on April 2, 2019 by Kennedy Mwangi
#
module load htseq/0.11.2

mkdir -p ~/tbrucei_gcn/results/HTSeq_count_results

GFF_FILE=$1

for sam_file in $(ls ~/tbrucei_gcn/data/raw_data/*.sam); do
    sam_file_name=$(echo $sam_file | cut -f1 -d '.')
    
        htseq-count \
            --format=sam \
            --stranded=yes \ #TODO: check whether stranded or not.
            --type=exon \ #feature type
            --idattr=Parent \
            $sam_file \
            $GFF_FILE > ${sam_file_name}.counts.txt
done
