# for changing transcript ids to corresponding gene ids
gtf_file <- import("../data/TriTrypDB-43_TbruceiTREU927.gtf")

gene_and_transcript_id <- mcols(gtf_file)[,c("gene_id","transcript_id")]

gene_and_transcript_id <- unique(gene_and_transcript_id)


# Remove extra column from sample SRR965341 added while reading data into R
data.all["Var.2"] <- NULL

# Create a DGEList object
counts <- DGEList(data.all, group = sample.metadata$Tissue)

# check the number of genes with no expression in all samples
table(rowSums(counts$counts==0)==15)
#FALSE  TRUE 
# 9184   792

# Filtering non-expressed and lowly-expressed genes.
keep.exprs <- filterByExpr(counts, group=sample.metadata$Sample_Name)
filtered.counts <- counts[keep.exprs,, keep.lib.sizes=FALSE]

# replace transcript ids with gene ids as rownames
filtered.counts.tmp <- tibble::rownames_to_column(as.data.frame(filtered.counts$counts),
                                                  "transcript_id")

filtered.counts.tmp$gene_id <- gene_and_transcript_id$gene_id[match(filtered.counts.tmp$transcript_id,
                                                                    gene_and_transcript_id$transcript_id)]

filtered.counts.tmp <- as.data.frame(filtered.counts.tmp) %>% remove_rownames %>% 
  column_to_rownames(var = "gene_id")

filtered.counts.tmp$transcript_id <- NULL

filtered.counts$counts <- as.matrix(filtered.counts.tmp)

# obtain logCPM unnormalized for plotting purposes.
# Here, the norm.factors value is 1 for all samples
logcpm.unnorm.counts <- cpm(filtered.counts, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

# Normalize for composition bias using TMM
filtered.counts <- calcNormFactors(filtered.counts, method = 'TMM')

# Convert counts per million per gene to log counts per million for further downstream analysis.
logcpm.norm.counts <- cpm(filtered.counts, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

# use ComBat to remove batch effects
modcombat <- model.matrix(~Tissue, data=sample.metadata)
logcpm.norm.counts.combat <- ComBat(dat=logcpm.norm.counts, batch = sample.metadata$Batch, 
                                    mod = modcombat)




#Various plots are made for the samples before and after normalization.

################################
#**Samples heatmap**
###############################

sample_category <- nlevels(sample.metadata$Sample_Name)
colour.palette <- colorRampPalette(brewer.pal(sample_category, "Set2"))(sample_category)
sample.colours <- colour.palette[as.integer(sample.metadata$Sample_Name)]

# Unnormalized sample heatmap
png(filename = "../figures/figure01_unnormalized-sample-heatmap.png", res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.unnorm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample")
dev.off()

# Normalized sample heatmap
png(filename = "../figures/figure01_normalized-sample-heatmap.png", res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample")

dev.off()

# plot to check removal of batch effects by ComBat
png(filename = "../figures/batch_effect_removal_analysis/figure01_combat-norm-sample-heatmap.png", res =1200, 
    type = "cairo", 
    units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts.combat), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample")

dev.off()



#################################
# **Samples density plot**
################################


# Checking further using sample density plot

# raw data
log.counts <- log2(counts$counts + 1)
png("../figures/figure02_raw-sample-density.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 6, pointsize = 10)
x <- melt(as.matrix(log.counts))

colnames(x) <- c('gene_id', 'sample', 'log')
ggplot(x, aes(x=log, color=sample)) + geom_density()
dev.off()

# filtered and unnormalized sample data
png("../figures/figure02_unnormalized-sample-density.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
x <- melt(as.matrix(logcpm.unnorm.counts))

colnames(x) <- c('gene_id', 'sample', 'logcpm')
ggplot(x, aes(x=logcpm, color=sample)) + geom_density()

dev.off()

# filtered and normalized sample data
png("../figures/figure02_normalized-sample-density.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
x <- melt(as.matrix(logcpm.norm.counts))

colnames(x) <- c('gene_id', 'sample', 'logcpm')
ggplot(x, aes(x=logcpm, color=sample)) + geom_density()

dev.off()

###########################################
# **Principal component analysis**
###########################################
# PCA

# raw samples PCA
pca.log.counts <- prcomp(t(log.counts)) # raw data (unnormalized and unfiltered)
png(filename = "../figures/figure03_raw-samples-PCA.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()

# unnormalized samples PCA
pca.log.counts <- prcomp(t(logcpm.unnorm.counts))  #unnormalized & filtered
png(filename = "../figures/figure03_unnormalized-sample-PCA.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()


# normalized samples PCA
pca.log.counts <- prcomp(t(logcpm.norm.counts))  #normalized & filtered
png(filename = "../figures/figure03_normalized-sample-PCA.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()

# plot to check removal of batch effects by ComBat
pca.log.counts.combat <- prcomp(t(logcpm.norm.counts.combat))  #normalized & filtered
png(filename = "../figures/batch_effect_removal_analysis/figure03_combat-normalized-sample-PCA.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts.combat,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()

#######################
#  **Boxplot**
#######################

# raw sample boxplot
png(filename = "../figures/figure04_raw-sample-boxplot.png", res =1200, type = "cairo", units = 'in',
    width = 4, height = 4, pointsize = 6)
y <- melt(as.matrix(log.counts))

colnames(y) <- c('gene_id', 'sample', 'log')
ggplot(y, aes(x=sample, y=log)) + geom_boxplot() + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()

# unnormalized sample boxplot
png(filename = "../figures/figure04_unnormalized-sample-boxplot.png", res =1200, type = "cairo", units = 'in',
    width = 4, height = 4, pointsize = 6)
y <- melt(as.matrix(logcpm.unnorm.counts))

colnames(y) <- c('gene_id', 'sample', 'logcpm')
ggplot(y, aes(x=sample, y=logcpm)) + geom_boxplot() + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()

# normalized sample boxplot
png(filename = "../figures/figure04_normalized-sample-boxplot.png", res =1200, type = "cairo", units = 'in',
    width = 4, height = 4, pointsize = 6)
y <- melt(as.matrix(logcpm.norm.counts))

colnames(y) <- c('gene_id', 'sample', 'logcpm')
ggplot(y, aes(x=sample, y=logcpm)) + geom_boxplot() + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()


