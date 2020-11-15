################################
#**Samples heatmap**
###############################

source(here::here("scripts","analysis","libraries.R"))
logcpm.unnorm.counts <- readRDS(here::here("data","intermediate","logcpm.unnorm.counts.RDS"))
logcpm.norm.counts <- readRDS(here::here("data","intermediate","logcpm.norm.counts.RDS"))
logcpm.norm.counts.combat <- readRDS(here::here("data","intermediate","logcpm.norm.counts.combat.RDS"))
samples.metadata.clean <- readRDS(file = here::here("data","raw","samples.metadata.clean.RDS"))

samples.metadata.clean$Sample_Name <- as.factor(samples.metadata.clean$Sample_Name)

sample_category <- nlevels(samples.metadata.clean$Sample_Name)
colour.palette <- colorRampPalette(brewer.pal(sample_category, "Set2"))(sample_category)
sample.colours <- colour.palette[as.integer(samples.metadata.clean$Sample_Name)]

# Raw sample heatmap
png(filename = here::here("results","figures","raw-sample-heatmap.png"), res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.unnorm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = samples.metadata.clean$Sample_Name,
          labCol = samples.metadata.clean$Sample_Name)
dev.off()

# Normalized sample heatmap
png(filename = here::here("results","figures","norm-sample-heatmap.png"), res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = samples.metadata.clean$Sample_Name,
          labCol = samples.metadata.clean$Sample_Name)

dev.off()

# plot to check removal of batch effects by ComBat
png(filename = here::here("results","figures","combat-norm-sample-heatmap.png"), res =1200, 
    type = "cairo", 
    units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts.combat), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = samples.metadata.clean$Sample_Name,
          labCol = samples.metadata.clean$Sample_Name)

dev.off()
