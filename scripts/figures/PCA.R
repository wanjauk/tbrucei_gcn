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