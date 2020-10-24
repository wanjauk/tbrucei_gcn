# Parameters
bam_file <- "../data/processed_data/bru-mor_bam/SRR039951.dupMarked.bam"
gtf_file <- "../data/brucei-morsitans/brucei-morsitans_annotations.gtf"
stranded <- 0
paired <- FALSE
threads <- 12

# Duplication rate ananlysis
dm <- analyzeDuprates(bam_file, gtf_file, stranded, paired, threads)

#Plots
png(filename = "../figures/duplication_rate/SRR039951.png")
duprateExpDensPlot(DupMat=dm)
title("SRR039951")
dev.off()

# Boxplot
duprateExpBoxplot(DupMat=dm)


# Exclusion of the following samples was done after analysis showed they
# failed quality control.

# remove 3 samples that have technical duplicates (SRR039951, SRR039937, SRR039938)
sample.metadata <- sample.metadata[-15,]
sample.metadata <- sample.metadata[-9,]
sample.metadata <- sample.metadata[-8,]