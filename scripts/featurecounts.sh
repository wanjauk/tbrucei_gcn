#!/bin/bash
#
#Script to counts reads that aligned to T. brucei genome using 
#featureCounts program from Subread package.
#
#Created April 10, 2019 by Kennedy Mwangi
#
#NOTE: featureCounts program needs to be installed on the HPC and
#the path to the program added to PATH
#
SAM_FILES=~/tbrucei_gcn/data/processed_data/*.sam

ANNOTATION_FILE=$1

for sam_file in $SAM_FILES; do
    sam_file_name=$(echo $sam_file | cut -f1 -d '.')

    /home/wanjau/subread-1.6.4-Linux-x86_64/bin/featureCounts \
        -a $ANNOTATION_FILE \
        -t exon \
        -g Parent \
        -T 1 \
        -o ${sam_file_name}.txt \
        $sam_file
done
