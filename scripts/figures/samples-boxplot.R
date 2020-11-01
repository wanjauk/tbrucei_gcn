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


