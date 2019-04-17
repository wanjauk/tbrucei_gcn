#!/bin/bash
#
#Script to counts reads that aligned to T. brucei genome using 
#featureCounts program from Subread package.
#
#Create output directory if it doesn't exist
mkdir -p ../results/featureCounts_results

SAM_FILES=../data/processed_data/*.sam

ANNOTATION_FILE=$1

for sam_file in $SAM_FILES; do
    sam_file_name=$(echo $sam_file | cut -f1 -d '.')

    featureCounts \
        -a $ANNOTATION_FILE \
        -t exon \
        -g Parent \
        -T 8 \
        -o ${sam_file_name}.txt \
        $sam_file
done

#move output files to results directory
#TODO: integrate the path in commandline option -o.
mv ../data/processed_data/SRR*.txt ../results/featureCounts_results/
