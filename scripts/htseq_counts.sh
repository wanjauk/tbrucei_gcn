#!/bin/bash
#
#Script to counts the number of reads aligned to T. brucei genome using HTSeq.
#Resource: HTSeq documentation https://htseq.readthedocs.io/en/latest/count.html
#Created on April 2, 2019 by Kennedy Mwangi
#
module load htseq/0.11.2

#create output directory if it doesn't exist
mkdir -p ../results/HTSeq_count_results

GFF_FILE=$1

for sam_file in ../data/processed_data/*.sam; do
    sam_file_name=$(echo $sam_file | cut -f1 -d '.')
    
        python /opt/apps/htseq/0.11.2/bin/htseq-count \
            -f sam \
            -s no \
            -t exon \
            -i Parent \
            $sam_file \
            $GFF_FILE \
            > ../results/HTSeq_count_results/${sam_file_name}.counts.txt
done
