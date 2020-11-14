# #!/bin/bash
# script to exclude gene features from reads counts 

# USAGE:
# ./exclude-features.sh \
# ../../data/intermediate/tbrucei_read_counts/savage \
# ../../data/intermediate/excluded_features.txt \
# ../../data/intermediate/tbrucei_read_counts_mRNA-only/savage \
# ../../data/intermediate/tbrucei_read_counts_mRNA-only/savage-telleria_tbrucei_read_counts_mRNA-only

# ./exclude-features.sh \
# ../../data/intermediate/tbrucei_read_counts/telleria \
# ../../data/intermediate/excluded_features.txt \
# ../../data/intermediate/tbrucei_read_counts_mRNA-only/telleria \
# ../../data/intermediate/tbrucei_read_counts_mRNA-only/savage-telleria_tbrucei_read_counts_mRNA-only

# create directory for the filtered counts
mkdir -p ../../data/intermediate/tbrucei_read_counts_mRNA-only/savage
mkdir -p ../../data/intermediate/tbrucei_read_counts_mRNA-only/telleria

# reads counts directory
READ_COUNTS_DIR=$1

# excludes features file
EXCLUDED_FEAT=$2

# mRNA only output directory
COUNTS_OUT=$3

#combined directory
COMBINED_DIR=$4

for file in ${READ_COUNTS_DIR}/*.counts.txt; do
  counts_file=$(basename "$file" .counts.txt)
		grep -v -f ${EXCLUDED_FEAT} ${file} > \
			${COUNTS_OUT}/${counts_file}.counts.txt
done


# copy savage and telleria read counts in a combined directory
# exclude samples which failed quality control (SRR039951, SRR039937, SRR039938)
mkdir -p $COMBINED_DIR

cp $COUNTS_OUT/*.counts.txt $COMBINED_DIR/
rm $COMBINED_DIR/SRR039951.counts.txt $COMBINED_DIR/SRR039937.counts.txt $COMBINED_DIR/SRR039938.counts.txt




