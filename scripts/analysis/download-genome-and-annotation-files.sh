#Downloading T. brucei genome

wget https://tritrypdb.org/common/downloads/release-43/TbruceiTREU927/fasta/data/TriTrypDB-43_TbruceiTREU927_Genome.fasta \
-P ../../data/scratch/tbrucei/

#Downloading the GFF file
wget https://tritrypdb.org/common/downloads/release-43/TbruceiTREU927/gff/data/TriTrypDB-43_TbruceiTREU927.gff \
-P ../../data/scratch/tbrucei/

# convert the tbrucei gene annotation from GFF format to GTF (required by some downstream tools)
# uses gffread from cufflinks
gffread ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.gff \
-T -o ../../data/scratch/tbrucei/TriTrypDB-43_TbruceiTREU927.gtf

# Download T. brucei annotated transcripts (for use in UTR motif discovery) 
wget https://tritrypdb.org/common/downloads/release-43/TbruceiTREU927/fasta/data/TriTrypDB-43_TbruceiTREU927_AnnotatedTranscripts.fasta \
-P ../../data/scratch/tbrucei/

# Downloading Glossina genome --Moved to new loaction after VEuPathDB creation
# wget https://www.vectorbase.org/download/glossina-morsitans-yalescaffoldsgmory1fagz \
# -P ../../data/scratch/glossina/

# Downloading GTF file --Moved to new loaction after VEuPathDB creation
# wget https://www.vectorbase.org/download/glossina-morsitans-yalebasefeaturesgmory19gtfgz \
# -P ../../data/scratch/glossina/

# Downloading Glossina genome
wget https://vectorbase.org/common/downloads/Pre-VEuPathDB%20VectorBase%20files/Glossina-morsitans-Yale_SCAFFOLDS_GmorY1.fa.gz \
-P ../../data/scratch/glossina/

# Downloading GTF file
wget https://vectorbase.org/common/downloads/Pre-VEuPathDB%20VectorBase%20files/Glossina-morsitans-Yale_BASEFEATURES_GmorY1.9.gtf.gz \
-P ../../data/scratch/glossina/

# unzip Glossina genome file
gunzip ../../data/scratch/glossina/Glossina-morsitans-Yale_SCAFFOLDS_GmorY1.fa.gz

# unzip Glossina annotation file file
gunzip ../../data/scratch/glossina/Glossina-morsitans-Yale_BASEFEATURES_GmorY1.9.gtf.gz






