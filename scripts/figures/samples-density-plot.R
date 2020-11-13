################################
# sample density plot
################################

source(here::here("scripts","analysis","libraries.R"))
logcpm.unnorm.counts <- readRDS(here::here("data","intermediate","logcpm.unnorm.counts.RDS"))
logcpm.norm.counts <- readRDS(here::here("data","intermediate","logcpm.norm.counts.RDS"))
reads_count <- readRDS(file = here::here("data", "intermediate", "tbrucei_reads_count.RDS"))

# raw data
log.counts <- log2(reads_count + 1)
png(here::here("results","figures","raw-sample-density.png"), res =1200, type = "cairo", units = 'in',
    width = 6, height = 6, pointsize = 10)
x <- melt(as.matrix(log.counts))

colnames(x) <- c('gene_id', 'sample', 'log')
ggplot(x, aes(x=log, color=sample)) + geom_density()
dev.off()

# # filtered and unnormalized sample data
# png(here::here("results","figures","unnormalized-sample-density.png"), res =1200, type = "cairo", units = 'in',
#     width = 6, height = 4, pointsize = 10)
# x <- melt(as.matrix(logcpm.unnorm.counts))
# 
# colnames(x) <- c('gene_id', 'sample', 'logcpm')
# ggplot(x, aes(x=logcpm, color=sample)) + geom_density()
# 
# dev.off()

# filtered and normalized sample data
png(here::here("results","figures","normalized-sample-density.png"), res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
x <- melt(as.matrix(logcpm.norm.counts))

colnames(x) <- c('gene_id', 'sample', 'logcpm')
ggplot(x, aes(x=logcpm, color=sample)) + geom_density()

dev.off()
