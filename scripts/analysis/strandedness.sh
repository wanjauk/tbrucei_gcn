# convert GTF genome annotation to BED format using a custom script from:
#https://github.com/ExpressionAnalysis/ea-utils/tree/master/clipper

../scripts/gtf2bed.pl ../data/tbrucei_genome_annotations_GTF/TriTrypDB-43_TbruceiTREU927.gtf > \
../data/tbrucei_genome_annotations_GTF/TriTrypDB-43_TbruceiTREU927.bed

infer_experiment.py -i ../data/processed_data/bru-mor_bam/SRR039381.bam \
-r ../data/tbrucei_genome_annotations_GTF/TriTrypDB-43_TbruceiTREU927.bed

# output
#This is SingleEnd Data
#Fraction of reads failed to determine: 0.0008
#Fraction of reads explained by "++,--": 0.3832
#Fraction of reads explained by "+-,-+": 0.6161