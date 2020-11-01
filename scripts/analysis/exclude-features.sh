# #!/bin/bash
# script to exclude gene features from reads counts 

# USAGE:
# ./exclude-features.sh \
# ../../data/intermediate/tbrucei_read_counts/savage \
# ../../data/intermediate/excluded_features.txt \
# ../../data/intermediate/tbrucei_read_counts_mRNA-only/savage

# create directory for the filtered counts
mkdir -p ../../data/intermediate/tbrucei_read_counts_mRNA-only/savage

# reads counts directory
READ_COUNTS_DIR=$1

# excludes features file
EXCLUDED_FEAT=$2

# mRNA only output directory
COUNTS_OUT=$3

for file in ${READ_COUNTS_DIR}/*.counts.txt; do
  counts_file=$(basename "$file" .counts.txt)
		grep -v -f ${EXCLUDED_FEAT} ${file} > \
			${COUNTS_OUT}/${counts_file}.counts.txt
done
