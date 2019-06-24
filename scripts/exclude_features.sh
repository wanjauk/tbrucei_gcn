#!/bin/bash
# script to exclude gene features from reads counts 

# create directory for the filtered counts
mkdir -p ../results/brucei_HTSeq_count_results_mRNA

for file in ../results/brucei_HTSeq_count_results/*.counts.txt; do
  counts_file=$(basename "$file" .counts.txt)
		grep -v -f ../results/excluded_features.txt ${file} > \
			../results/brucei_HTSeq_count_results_mRNA/"${counts_file}".counts.txt
done
