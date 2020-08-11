# July 2, 2020
#########################
## Sample heatmap
#########################
sample.metadata$Sample_Name <- as.factor(sample.metadata$Sample_Name)

sample_category <- nlevels(sample.metadata$Sample_Name)
colour.palette <- colorRampPalette(brewer.pal(sample_category, "Set2"))(sample_category)
sample.colours <- colour.palette[as.integer(sample.metadata$Sample_Name)]

# Unnormalized sample heatmap
png(filename = "../figures/unnorm_sample_heatmap_july_2020.png", res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.unnorm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = sample.metadata$Sample_Name,
          labCol = sample.metadata$Sample_Name)
dev.off()

# Normalized sample heatmap
png(filename = "../figures/norm_sample_heatmap_july_2020.png", res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = sample.metadata$Sample_Name,
          labCol = sample.metadata$Sample_Name)

dev.off()

# plot to check removal of batch effects by ComBat
png(filename = "../figures/batch_effect_removal_analysis/combat_norm_sample_heatmap_july_2020.png", res =1200, 
    type = "cairo", 
    units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts.combat), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = sample.metadata$Sample_Name,
          labCol = sample.metadata$Sample_Name)

dev.off()


#############################################
# PCA
###############################################

# raw samples PCA
# pca.log.counts <- prcomp(t(log.counts)) # raw data (unnormalized and unfiltered)
# 
# png(filename = "../figures/raw_samples_PCA.png", res =1200, type = "cairo", units = 'in',
#     width = 6, height = 4, pointsize = 10)
# autoplot(pca.log.counts,
#          data = sample.metadata,
#          colour="Sample_Name",
#          size=3)
# dev.off()

# unnormalized samples PCA
pca.log.counts <- prcomp(t(logcpm.unnorm.counts))  #unnormalized & filtered
png(filename = "../figures/unnorm_sample_PCA_july_2020.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()


# normalized samples PCA
pca.log.counts <- prcomp(t(logcpm.norm.counts))  #normalized & filtered
png(filename = "../figures/norm_sample_PCA_july_2020.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()

# plot to check removal of batch effects by ComBat
pca.log.counts.combat <- prcomp(t(logcpm.norm.counts.combat))  #normalized & filtered
png(filename = "../figures/batch_effect_removal_analysis/combat_norm_sample_PCA_july_2020.png", res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts.combat,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()