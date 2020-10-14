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
    bam_file_name=$(basename "$bam_file" .bam)
    
        python /opt/apps/htseq/0.11.2/bin/htseq-count \
            -f bam \
            -s no \
            -t exon \
            -i Parent \
            $bam_file \
            $GFF_FILE \
            > ../results/brucei_HTSeq_count_results/${bam_file_name}.counts.txt
done



## Generating MultiQC report

#MultiQC aggregates results from FASTQC, HISAT2 and HTSeq analysis into an HTML formatted single report for better visualization.
#change directory to results
cd ../results

#Run multiqc
multiqc .

# create a directory for the multiQC report and move the output there.
mkdir -p brucei_multiqc_report
mv multiqc* brucei_multiqc_report/