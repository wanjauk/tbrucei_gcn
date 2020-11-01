# Load SRA metadata from Savage et al 2016 study
SRA.metadata <- read.table(here::here("data", "raw", "savage","SraRunTable.metadata.txt"), 
                           header = TRUE, sep = "\t")

# obtain sample metadata to be used later in analysis in R.
matches <- c("Run","Library_Name","Sample_Name")
samples.metadata <- SRA.metadata[grepl(paste(matches, collapse="|"), names(SRA.metadata))]

# create a grouping factor that will place each sample in the one of three tissues i.e.
# midgut (MG), proventriculus(PV) and salivary glands (SG)
tissue <- factor(c("MG", "MG", "MG", "MG", "MG", "PV", "PV", "SG", "SG", "SG", "SG", 
                   "MG", "MG", "PV", "SG", "SG", "PV"))

# append factor to samples.metadata to group samples
samples.metadata["Tissue"] <- tissue

# Add sample from Telleria et al 2014 study (SRR965341) to sample metadata.
samples.metadata$Run <- as.character(samples.metadata$Run)
samples.metadata <- rbind(samples.metadata, "18" = c("SA2", "SRR965341", "SA2", "SG"))

# include batch information for the 18 samples
batch <- factor(c(1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 2, 1, 2, 2, 2))
samples.metadata["Batch"] <- batch

# save the metadata to an R object
saveRDS(samples.metadata, here::here("data", "raw", "samples.metadata.RDS"))
        
        