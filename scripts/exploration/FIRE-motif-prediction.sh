
#Prepare input sequences for FIRE by obtaining module genes sequences from *T. brucei*'s entire annotated transcripts. Use `seqtk` package.


# remove unrequired info from header and retain gene id.
cut -f2 -d "|" ../data/tbrucei_annotated_transcripts/TriTrypDB-43_TbruceiTREU927_AnnotatedTranscripts.fasta | cut -f2 -d " " | cut -f1 -d " " | cut -f2 -d "=" > ../results/tbrucei_annotated_transcript_sequences_gene_as_header.fasta

# add ">" symbol to header
awk '{ if ($0 ~ /Tb/) { printf ">"; } print $0; }' ../results/tbrucei_annotated_transcript_sequences_gene_as_header.fasta > ../results/tbrucei_annotated_transcript_sequences.fasta

# obtain module gene sequences
seqtk subseq ../results/tbrucei_annotated_transcript_sequences.fasta \
../results/tbrucei_module_genes.txt > \
../results/tbrucei_module_genes_sequences_FIRE_input.fasta

seqkit grep -n -f ../results/tbrucei_module_genes_sorted.txt ../results/tbrucei_annotated_transcript_sequences.fasta > ../results/tbrucei_module_genes_sequences_FIRE_input.fasta1

# FIRE commands - FIRE should be run in the directory where all scripts were installed. Therefore,
# the input files should be copied to FIRE's directory (FIRE-x.x/)
#
#check whether input files are OK
perl TOOLS/FIRE_analyse_input_files.pl \
-fastafile tbrucei_module_genes_sequences_without_delimiter.fasta \
-expfile tbrucei_FIRE_expression_clusters.txt

# remove sequences >= 10,000 base pairs as they make FIRE crash
# code - https://www.biostars.org/p/62678/
 cat ../results/tbrucei_module_genes_sequences_without_delimiter.fasta | \
 awk '{y= i++ % 2 ; L[y]=$0; if(y==1 && length(L[1])<=9999) {printf("%s\n%s\n",L[0],L[1]);}}' > \
../results/tbrucei_module_genes_sequences_without_delimiter_filtered.fasta

# recheck input files
perl TOOLS/FIRE_analyse_input_files.pl \
-fastafile tbrucei_module_genes_sequences_without_delimiter_filtered.fasta \
-expfile tbrucei_FIRE_expression_clusters.txt

# output
#Checking the expression file
#Expression file is OK.
#
#Checking the fasta file
#Fasta file is OK (min sequence length = 98, max length = 9917).
#
#Found fasta sequence for 7272 / 7390 identifiers in expression file.

# Analyze sequences for motifs
 perl fire.pl --expfiles=tbrucei_FIRE_expression_clusters.txt --exptype=discrete \
 --fastafile_dna=tbrucei_module_genes_sequences_without_delimiter_filtered.fasta  --nodups=1

# create HTML output of the results
perl MORESCRIPTS/makeresultindex.pl tbrucei_FIRE_expression_clusters.txt "T. brucei cluster motifs"

