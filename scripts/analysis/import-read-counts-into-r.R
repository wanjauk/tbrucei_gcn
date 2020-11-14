#!/usr/bin/Rscript

# Take 'all' htseq-count results and melt them in to one big dataframe

#Adapted from: https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4

# required packages
library(tibble)

# # where are we?
cntdir <- here::here("data", "intermediate", "tbrucei_read_counts_mRNA-only",
                    "savage-telleria_tbrucei_read_counts_mRNA-only")
pat <- ".counts.txt"
hisat2.all <- list.files(path = cntdir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# we choose the 'all' series
myfiles <- hisat2.all
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*).counts.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
tbrucei_reads_count_raw <- DT[[myfiles[1]]]

# inspect
#head(tbrucei_reads_count_raw)

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(tbrucei_reads_count_raw, y, by = c("ID"))
  tbrucei_reads_count_raw <- z
}

# ID column becomes rownames
rownames(tbrucei_reads_count_raw) <- tbrucei_reads_count_raw$ID
tbrucei_reads_count_raw <- tbrucei_reads_count_raw[,-1]

## add total counts per sample
tbrucei_reads_count_raw <- rbind(tbrucei_reads_count_raw, 
                                 tot.counts=colSums(tbrucei_reads_count_raw))

# inspect and look at the top row names!
#head(tbrucei_reads_count_raw)

#tail(tbrucei_reads_count_raw)

####################################
# take summary rows to a new table
# ( not starting with Tb and tmp with invert=TRUE )

# transpose table for readability
reads_count_summary <- tbrucei_reads_count_raw[grep("^Tb|^tmp", rownames(tbrucei_reads_count_raw), 
                                                    perl=TRUE, invert=TRUE), ]

# review
#reads_count_summary

# transpose table
t(reads_count_summary)

# write summary to file
write.csv(reads_count_summary, 
          file = here::here("data", "intermediate","tbrucei_reads_count_summary.csv"), 
          row.names = FALSE)

####################################
# take all data rows to a new table
reads_count <- tbrucei_reads_count_raw[grep("^Tb|^tmp", rownames(tbrucei_reads_count_raw), perl=TRUE, invert=FALSE), ]

# inspect final merged table
#head(reads_count, 3)

# write data to files
saveRDS(reads_count, file = here::here("data", "intermediate", "tbrucei_reads_count.RDS"))

reads_count <- rownames_to_column(reads_count,"transcript_id")
write.csv(reads_count, 
          file = here::here("data", "intermediate", "tbrucei_reads_count.csv"), 
          row.names = FALSE)

# cleanup intermediate objects
rm(y, z, i, DT, tbrucei_reads_count_raw)
