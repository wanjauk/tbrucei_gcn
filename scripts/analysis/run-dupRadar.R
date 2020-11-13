load(file = here::here("data","raw","sample.metdata.final.RData"))
source(here::here("scripts","analysis","libraries.R"))
source(here::here("scripts","analysis","settings.R"))

# Parameters
gtf_file <- here::here("data","scratch","concatenated_genomes","brucei-morsitans_annotations.gtf")
stranded <- 0
paired <- FALSE
threads <- 20

# path to bam files data
bam_file <- list.files(here::here("data","scratch","mark-duplicates-output"),
                       pattern = ".bam",
                       full.names = TRUE)

for (file in bam_file){
  
  #get the file name
  file_name <- gsub(pattern = "\\.bam$", "", file)
  
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



# Exclusion of the following samples was done after analysis showed they
# failed quality control.

# remove 3 samples that have technical duplicates (SRR039951, SRR039937, SRR039938)
# sample.metadata <- sample.metadata[-15,]
# sample.metadata <- sample.metadata[-9,]
# sample.metadata <- sample.metadata[-8,]
# 
saveRDS(sample.metadata, here::here("data","raw","sample.metadata.RDS"))

