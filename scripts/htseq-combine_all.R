#!/usr/bin/Rscript

# Take 'all' htseq-count results and melt them in to one big dataframe

#Adapted from: https://wiki.bits.vib.be/index.php/NGS_RNASeq_DE_Exercise.4

# where are we?
basedir <- "../results"
setwd(basedir)

cntdir <- paste(basedir, "brucei_HTSeq_count_results", sep="/")
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
data <- DT[[myfiles[1]]]

# inspect
head(data)

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))

# inspect and look at the top row names!
head(data)

tail(data)

####################################
# take summary rows to a new table
# ( not starting with Tb and tmp with invert=TRUE )

# transpose table for readability
data.all.summary <- data[grep("^Tb|^tmp", rownames(data), perl=TRUE, invert=TRUE), ]

# review
data.all.summary

# transpose table
t(data.all.summary)

# write summary to file
write.csv(data.all.summary, file = "brucei_htseq_counts_all-summary.csv")

####################################
# take all data rows to a new table

data.all <- data[grep("^Tb|^tmp", rownames(data), perl=TRUE, invert=FALSE), ]

# inspect final merged table
head(data.all, 3)

# write data to file
write.table(data.all, file = "brucei_htseq_counts_all.txt", quote = FALSE, sep = "\t")

# cleanup intermediate objects
rm(y, z, i, DT)