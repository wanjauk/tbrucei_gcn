################################
#**Samples heatmap**
###############################

source(here::here("scripts","analysis","libraries.R"))
logcpm.unnorm.counts <- readRDS(here::here("data","intermediate","logcpm.unnorm.counts.RDS"))
logcpm.norm.counts <- readRDS(here::here("data","intermediate","logcpm.norm.counts.RDS"))
logcpm.norm.counts.combat <- readRDS(here::here("data","intermediate","logcpm.norm.counts.combat.RDS"))
load(file = here::here("data","raw","sample.metdata.final.RData"))

sample.metadata$Sample_Name <- as.factor(sample.metadata$Sample_Name)

sample_category <- nlevels(sample.metadata$Sample_Name)
colour.palette <- colorRampPalette(brewer.pal(sample_category, "Set2"))(sample_category)
sample.colours <- colour.palette[as.integer(sample.metadata$Sample_Name)]

# Unnormalized sample heatmap
png(filename = here::here("results","figures","unnorm-sample-heatmap"), res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.unnorm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = sample.metadata$Sample_Name,
          labCol = sample.metadata$Sample_Name)
dev.off()

# Normalized sample heatmap
png(filename = here::here("results","figures","norm-sample-heatmap.png"), res =1200, type = "cairo", units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = sample.metadata$Sample_Name,
          labCol = sample.metadata$Sample_Name)

dev.off()

# plot to check removal of batch effects by ComBat
png(filename = here::here("results","figures","combat-norm-sample-heatmap.png"), res =1200, 
    type = "cairo", 
    units = 'in',
    width = 5, height = 4, pointsize = 10)
heatmap.2(cor(logcpm.norm.counts.combat), RowSideColors=sample.colours, trace='none', 
          main='Sample correlations', margins = c(8, 8), xlab="Sample", ylab="Sample",
          labRow = sample.metadata$Sample_Name,
          labCol = sample.metadata$Sample_Name)

dev.off()
