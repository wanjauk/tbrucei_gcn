source(here::here("scripts","analysis","libraries.R"))
logcpm.unnorm.counts <- readRDS(here::here("data","intermediate","logcpm.unnorm.counts.RDS"))
logcpm.norm.counts <- readRDS(here::here("data","intermediate","logcpm.norm.counts.RDS"))
reads_count <- readRDS(file = here::here("data", "intermediate", "tbrucei_reads_count.RDS"))

#######################
#  **Boxplot**
#######################

# raw sample boxplot
log.counts <- log2(reads_count + 1)
png(filename = here::here("results","figures","raw-sample-boxplot.png"), res =1200, type = "cairo", units = 'in',
    width = 4, height = 4, pointsize = 6)
y <- melt(as.matrix(log.counts))

colnames(y) <- c('gene_id', 'sample', 'log')
ggplot(y, aes(x=sample, y=log)) + geom_boxplot() + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()

# # unnormalized sample boxplot
# png(filename = here::here("results","figures","unnormalized-sample-boxplot.png"), res =1200, type = "cairo", units = 'in',
#     width = 4, height = 4, pointsize = 6)
# y <- melt(as.matrix(logcpm.unnorm.counts))
# 
# colnames(y) <- c('gene_id', 'sample', 'logcpm')
# ggplot(y, aes(x=sample, y=logcpm)) + geom_boxplot() + 
#   theme(axis.text.x  = element_text(angle=90, vjust=0.5))
# dev.off()

# normalized sample boxplot
png(filename = here::here("results","figures","normalized-sample-boxplot.png"), res =1200, type = "cairo", units = 'in',
    width = 4, height = 4, pointsize = 6)
y <- melt(as.matrix(logcpm.norm.counts))

colnames(y) <- c('gene_id', 'sample', 'logcpm')
ggplot(y, aes(x=sample, y=logcpm)) + geom_boxplot() + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()


