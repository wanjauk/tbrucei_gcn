#!/bin/bash
#
#Script to counts the number of reads aligned to T. brucei genome using HTSeq.
#Resource: HTSeq documentation https://htseq.readthedocs.io/en/latest/count.html
#
module load htseq/0.11.2

#create output directory if it doesn't exist
mkdir -p ../results/brucei_HTSeq_count_results

GFF_FILE=$1

for bam_file in ../data/processed_data/bru-mor_bam/*.bam; do
    bam_file_name=$(echo $bam_file | cut -f1 -d '.')
    
        python /opt/apps/htseq/0.11.2/bin/htseq-count \
            -f bam \
            -s no \
            -t exon \
            -i Parent \
            $bam_file \
            $GFF_FILE \
            > ../results/brucei_HTSeq_count_results/${bam_file_name}.counts.txt
done
