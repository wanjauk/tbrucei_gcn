###########################################
# **Principal component analysis**
###########################################

source(here::here("scripts","analysis","libraries.R"))
logcpm.unnorm.counts <- readRDS(here::here("data","intermediate","logcpm.unnorm.counts.RDS"))
logcpm.norm.counts <- readRDS(here::here("data","intermediate","logcpm.norm.counts.RDS"))
logcpm.norm.counts.combat <- readRDS(here::here("data","intermediate","logcpm.norm.counts.combat.RDS"))
load(file = here::here("data","raw","sample.metdata.final.RData"))
reads_count <- readRDS(file = here::here("data", "intermediate", "tbrucei_reads_count.RDS"))

# PCA

# raw samples PCA
log.counts <- log2(reads_count + 1)
pca.log.counts <- prcomp(t(log.counts)) # raw data (unnormalized and unfiltered)
png(filename = here::here("results","figures","raw-samples-PCA.png"), res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()

# # unnormalized samples PCA
# pca.log.counts <- prcomp(t(logcpm.unnorm.counts))  #unnormalized & filtered
# png(filename = here::here("results","figures","unnorm-sample-PCA.png"), res =1200, type = "cairo", units = 'in',
#     width = 6, height = 4, pointsize = 10)
# autoplot(pca.log.counts,
#          data = sample.metadata,
#          colour="Sample_Name",
#          size=3)
# dev.off()


# normalized samples PCA
pca.log.counts <- prcomp(t(logcpm.norm.counts))  #normalized & filtered
png(filename = here::here("results","figures","norm-sample-PCA.png"), res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()

# plot to check removal of batch effects by ComBat
pca.log.counts.combat <- prcomp(t(logcpm.norm.counts.combat))  #normalized & filtered
png(filename = here::here("results","figures","combat-norm-sample-PCA.png"), res =1200, type = "cairo", units = 'in',
    width = 6, height = 4, pointsize = 10)
autoplot(pca.log.counts.combat,
         data = sample.metadata,
         colour="Sample_Name",
         size=3)
dev.off()
