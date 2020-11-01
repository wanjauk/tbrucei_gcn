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
