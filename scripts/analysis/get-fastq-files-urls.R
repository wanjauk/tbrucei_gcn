#Obtain FASTQ urls to download the data used in this study from ENA database.

# ENA metadata
# code adopted from: 
# https://wiki.bits.vib.be/index.php/Download_read_information_and_FASTQ_data_from_the_SRA

# accession numbers from Savage et al and Telleria et al the studies
accessions <- c("SRP002243","SRR965341")

# create directories to store the files if they don't exist.
telleria_dir <- here::here("data", "raw", "telleria")
if (!dir.exists(telleria_dir)) {
  dir.create(telleria_dir, recursive=TRUE)
}

savage_dir <- here::here("data", "raw", "savage")
if (!dir.exists(savage_dir)) {
  dir.create(savage_dir, recursive=TRUE)
}

for (accession_num in accessions) {
  
  # these samples from Savage et al are single-end reads.
  if (accession_num == "SRP002243"){
    
    # construct the url to ENA database
    ena.url <- paste("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                     accession_num,
                     "&result=read_run",
                     "&fields=study_accession,sample_accession,",
                     "experiment_accession,run_accession,scientific_name,fastq_ftp,",
                     "submitted_ftp,sra_ftp",
                     "&download=true",
                     sep="")
    
    # Load the metadata from ENA
    ENA.metadata <- read.table(url(ena.url), header=TRUE, sep="\t")
    
    # get the fastq urls
    fastq.urls <- ENA.metadata[grepl("fastq_ftp", names(ENA.metadata))]
    
    # create a text file with urls to fastq files in ENA database
    write.table(fastq.urls, here::here(savage_dir, "savage.fastq.urls.txt"), 
                eol = "\n", 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
    
  } else {
    
    # The other sample (SRR965341) from Telleria et al study has paired-end reads.
    #
    # construct the url to ENA database
    ena.url <- paste("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",
                     accession_num,
                     "&result=read_run",
                     "&fields=study_accession,sample_accession,",
                     "experiment_accession,run_accession,scientific_name,fastq_ftp,",
                     "submitted_ftp,sra_ftp",
                     "&download=true",
                     sep="")
    
    # Load the metadata from ENA
    ENA.metadata <- read.table(url(ena.url), header=TRUE, sep="\t")
    
    # ensure that R1 and R2 fastq files are in separate rows
    ENA.metadata <- ENA.metadata %>% separate_rows(fastq_ftp, sep=";")
    
    # get the fastq urls
    fastq.urls <- ENA.metadata[grepl("fastq_ftp", names(ENA.metadata))]
    
    # create a text file with urls to fastq files in ENA database
    write.table(fastq.urls, here::here(telleria_dir, "telleria.fastq.urls.txt"), 
                eol = "\n", 
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }
  
}
