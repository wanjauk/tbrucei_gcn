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

