#Downloading T. brucei genome

wget https://tritrypdb.org/common/downloads/release-43/TbruceiTREU927/fasta/data/\
TriTrypDB-43_TbruceiTREU927_Genome.fasta \
-P ../data/tbrucei_genome/

#Downloading the GFF file
wget https://tritrypdb.org/common/downloads/release-43/TbruceiTREU927/gff/data/\
TriTrypDB-43_TbruceiTREU927.gff \
-P ../data/tbrucei_genome_annotations_GFF/

# convert the tbrucei gene annotation from GFF format to GTF (required by some downstream tools)
# uses gffread from cufflinks
mkdir -p ../data/tbrucei_genome_annotations_GTF

gffread ../data/tbrucei_genome_annotations_GFF/TriTrypDB-43_TbruceiTREU927.gff \
-T -o ../data/tbrucei_genome_annotations_GTF/TriTrypDB-43_TbruceiTREU927.gtf

# Download T. brucei annotated transcripts (for use in UTR motif discovery) 
wget https://tritrypdb.org/common/downloads/release-43/TbruceiTREU927/fasta/data/\
TriTrypDB-43_TbruceiTREU927_AnnotatedTranscripts.fasta \
-P ../data/tbrucei_annotated_transcripts/

# Downloading Glossina genome
wget https://www.vectorbase.org/download/glossina-morsitans-yalescaffoldsgmory1fagz \
-P ../data/glossina_genome_scaffolds/

# Downloading GTF file
wget https://www.vectorbase.org/download/glossina-morsitans-yalebasefeaturesgmory19gtfgz \
-P ../data/glossina_genome_annonations_GTF/
