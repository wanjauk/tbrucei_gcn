# convert GTF genome annotation to BED format using a custom script from:
#https://github.com/ExpressionAnalysis/ea-utils/tree/master/clipper/gtf2bed

# USAGE:
# # savage data
# ./check-reads-strandedness.sh \
# ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.gtf \
# ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.bed \
# ../../data/scratch/sam-to-bam-output/savage/SRR039381.sorted.bam

# # telleria data
# ./check-reads-strandedness.sh \
# ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.gtf \
# ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.bed \
# ../../data/scratch/sam-to-bam-output/telleria/SRR965341.sorted.bam

# GTF file
GTF_FILE=$1

# BED file directory
BED_FILE=$2

# BAM file
BAM_FILE=$3


# create a BED file from GTF
# ../utils/gtf2bed.pl $GTF_FILE > $BED_FILE

infer_experiment.py -i $BAM_FILE -r $BED_FILE

# output for savage's SRR039381.sorted.bam
# This is SingleEnd Data
# Fraction of reads failed to determine: 0.0000
# Fraction of reads explained by "++,--": 0.4058
# Fraction of reads explained by "+-,-+": 0.5942

# output for telleria's SRR965341.sorted.bam
# This is PairEnd Data
# Fraction of reads failed to determine: 0.0019
# Fraction of reads explained by "1++,1--,2+-,2-+": 0.4999
# Fraction of reads explained by "1+-,1-+,2++,2--": 0.4982
