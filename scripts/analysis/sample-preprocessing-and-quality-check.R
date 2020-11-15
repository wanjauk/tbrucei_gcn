# load required packages and data
source(here::here("scripts","analysis","libraries.R"))
samples.metadata.clean <- readRDS(here::here("data","raw","samples.metadata.clean.RDS"))
reads_count <- readRDS(file = here::here("data", "intermediate", "tbrucei_reads_count.RDS"))

# Create a DGEList object
counts <- DGEList(reads_count, group = samples.metadata.clean$Tissue)

# check the number of genes with no expression in all samples
table(rowSums(counts$counts==0)==15)
#FALSE  TRUE 
# 9184   792

# Filtering non-expressed and lowly-expressed genes.
keep.exprs <- filterByExpr(counts, group=samples.metadata.clean$Sample_Name)
filtered.counts <- counts[keep.exprs,, keep.lib.sizes=FALSE]

# # change transcript ids to corresponding gene ids -----------------------------------------
gtf_file <- import(here::here("data","scratch","tbrucei","TriTrypDB-43_TbruceiTREU927.gtf"))
gene_and_transcript_id <- mcols(gtf_file)[,c("gene_id","transcript_id")]
gene_and_transcript_id <- unique(gene_and_transcript_id)


# replace transcript ids with gene ids as rownames
filtered.counts.tmp <- tibble::rownames_to_column(as.data.frame(filtered.counts$counts),
                                                  "transcript_id")

filtered.counts.tmp$gene_id <- gene_and_transcript_id$gene_id[match(filtered.counts.tmp$transcript_id,
                                                                    gene_and_transcript_id$transcript_id)]

filtered.counts.tmp <- as.data.frame(filtered.counts.tmp) %>% remove_rownames %>%
  column_to_rownames(var = "gene_id")

filtered.counts.tmp$transcript_id <- NULL
filtered.counts$counts <- as.matrix(filtered.counts.tmp)

# obtain logCPM unnormalized for plotting purposes ------------------------------------------
# Here, the norm.factors value is 1 for all samples
logcpm.unnorm.counts <- cpm(filtered.counts, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

# Normalize for composition bias using TMM
filtered.counts <- calcNormFactors(filtered.counts, method = 'TMM')

# Convert counts per million per gene to log counts per million for further downstream analysis.
logcpm.norm.counts <- cpm(filtered.counts, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

# use ComBat to remove batch effects
modcombat <- model.matrix(~Tissue, data=samples.metadata.clean)
logcpm.norm.counts.combat <- ComBat(dat=logcpm.norm.counts, batch = samples.metadata.clean$Batch, 
                                    mod = modcombat)

# save outputs for later analysis
saveRDS(filtered.counts, here::here("data","intermediate","filtered.counts.RDS"))
saveRDS(logcpm.unnorm.counts, here::here("data","intermediate","logcpm.unnorm.counts.RDS"))
saveRDS(logcpm.norm.counts, here::here("data","intermediate","logcpm.norm.counts.RDS"))
saveRDS(logcpm.norm.counts.combat, here::here("data","intermediate","logcpm.norm.counts.combat.RDS"))
saveRDS(gene_and_transcript_id, here::here("data","intermediate","gene_and_transcript_id.RDS"))

# clean up the environment
rm(gtf_file, filtered.counts.tmp)

