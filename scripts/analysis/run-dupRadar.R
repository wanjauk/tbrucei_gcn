samples.metadata <- readRDS(file = here::here("data","raw","samples.metadata.RDS"))
source(here::here("scripts","analysis","libraries.R"))
source(here::here("scripts","analysis","settings.R"))

# Parameters
gtf_file <- here::here("data","scratch","concatenated_genomes",
                       "brucei-morsitans_annotations.gtf")
stranded <- 0
paired <- FALSE # for savage's data
threads <- 10

# path to bam files data
bam_file <- list.files(here::here("data","scratch","mark-duplicates-output","savage"),
                       pattern = ".bam",
                       full.names = TRUE)

for (file in bam_file){
      
  #get the file name
  file_name <- gsub(pattern = "\\.dupMarked.bam$", "", basename(file))

  # Duplication rate analysis
  dm <- analyzeDuprates(file, gtf_file, stranded, paired, threads)
  
  #Plots
  png(filename = here::here("results","figures","duplication_rate",paste0(file_name,".png")))
  duprateExpDensPlot(DupMat=dm)
  title(paste0(file_name))
  dev.off()
  
  # # Boxplot
  # duprateExpBoxplot(DupMat=dm)
  
}

#########################################
# Telleria paired-end data processing
#########################################

#Parameters
telleria_paired <- TRUE

# path to bam files data
telleria_bam_file <- list.files(here::here("data","scratch","mark-duplicates-output","telleria"),
                       pattern = ".bam",
                       full.names = TRUE)
#get the file name
telleria_file_name <- gsub(pattern = "\\.dupMarked.bam$", "", basename(telleria_bam_file))

# Duplication rate analysis
dm <- analyzeDuprates(telleria_bam_file, gtf_file, stranded, telleria_paired, threads)

#Plots
png(filename = here::here("results","figures","duplication_rate",paste0(telleria_file_name,".png")))
duprateExpDensPlot(DupMat=dm)
title(paste0(telleria_file_name))
dev.off()

# Exclusion of the following samples was done after analysis showed they
# failed quality control.

# remove 3 samples that have technical duplicates (SRR039951, SRR039937, SRR039938)
samples.metadata <- samples.metadata[-15,]
samples.metadata <- samples.metadata[-9,]
samples.metadata <- samples.metadata[-8,]

# save the samples metadata
saveRDS(samples.metadata, here::here("data","raw","samples.metadata.clean.RDS"))

