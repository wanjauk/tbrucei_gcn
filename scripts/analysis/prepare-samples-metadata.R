# SRA metadata
SRA.metadata <- read.table("../data/SraRunTable.metadata.txt", header = TRUE, sep = "\t")

# obtain sample metadata to be used later in analysis in R.
matches <- c("Run","Library_Name","Sample_Name")
sample.metadata <- SRA.metadata[grepl(paste(matches, collapse="|"), names(SRA.metadata))]

# create grouping factor that will place each sample in the one of three tissues i.e.
# midgut (MG), proventriculus(PV) and salivary glands (SG)
tissue <- factor(c("MG", "MG", "MG", "MG", "MG", "PV", "PV", "SG", "SG", "SG", "SG", 
                   "MG", "MG", "PV", "SG", "SG", "PV"))

# append factor to sample.metadata to group samples
sample.metadata["Tissue"] <- tissue

# The sample below was analysed separately as the reads are paired-end while
# the other samples are single-end.
#
# Add sample from Telleria et al 2014 study (SRR965341) to sample metadata.
sample.metadata$Run <- as.character(sample.metadata$Run)

sample.metadata <- rbind(sample.metadata, "18" = c("SA2", "SRR965341", "SA2", "SG"))

# include batch information for the samples
batch <- factor(c(1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2))

sample.metadata["Batch"] <- batch

