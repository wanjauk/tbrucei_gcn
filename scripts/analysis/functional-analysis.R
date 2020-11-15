source(here::here("scripts","analysis","libraries.R"))
source(here::here("scripts","analysis","settings.R"))
source(here::here("scripts","utils","enrichment_analysis.R"))
source(here::here("scripts","utils","annotations.R"))
source(here::here("scripts","utils","wgcna.R"))
source(here::here("scripts","utils","util.R"))
logcpm.norm.counts.combat <- readRDS(here::here("data","intermediate","logcpm.norm.counts.combat.RDS"))
gene_and_transcript_id <- readRDS(here::here("data","intermediate","gene_and_transcript_id.RDS"))
gene.tree <- readRDS(file = here::here("data","intermediate","gene.tree.RDS"))
all.modules <- readRDS(file = here::here("data","intermediate","all.modules.RDS"))
module.colours <- readRDS(file = here::here("data","intermediate","module.colours.RDS"))
module.hub.genes <- readRDS(here::here("data","intermediate","module.hub.genes.RDS"))

#################################################################
# load gene annotations from packages
#################################################################
orgdb <- get("Trypanosoma.brucei.TREU927")

# Fix AnnotationDbi namespace mess
assign('select', dplyr::select, envir=.GlobalEnv)
assign('get',    base::get, envir=.GlobalEnv)

gene_info <- load_parasite_annotations(orgdb, rownames(logcpm.norm.counts.combat),
                                       keytype="GID")

# add the transcript id column to the gene_info
gene_info$transcript_id <- gene_and_transcript_id$transcript_id[match(gene_info$gene_id,
                                                                    gene_and_transcript_id$gene_id)]

# Get transcript lengths (sum of all exon lengths for each gene)
txdb <- orgdb@txdbSlot
transcript_lengths <- transcriptLengths(txdb)
transcript_lengths <- transcript_lengths[transcript_lengths$tx_name %in%
                                           gene_info$transcript_id,]

# add the transcript lengths to gene_info
gene_info[match(transcript_lengths$tx_name, gene_info$transcript_id),
          'transcript_length'] <- transcript_lengths$tx_len
gene_info$transcript_length <- as.numeric(gene_info$transcript_length)

# A gene was excluded in the annotation database as a result of an orphan transcript.
#Add its placeholder.
## ---not run---
# gene_info <- rbind(gene_info, data.frame(gene_id='Tb927.4.4663',
#                                          chromosome='4',
#                                          description=NA, strand=NA, type=NA,
#                                          transcript_length=NA))


# Keep only the feature information remaining genes
gene_info <- gene_info[gene_info$gene_id %in% rownames(logcpm.norm.counts.combat),] #WGCNA input counts

# For now, just grab the description for the first transcript
#gene_info <- gene_info[!duplicated(gene_info$gene_id),]

# Gene IDs
gene_ids <- rownames(logcpm.norm.counts.combat)

# gene annotations preview
kable(head(gene_info), caption='Preview of gene annotations.')

#########################################################################################
# load GO terms associated with each parasite gene
#########################################################################################

# # load go terms from annotation package
# go_terms <- load_go_terms(orgdb, rownames(logcpm.norm.counts.combat),
#                           keytype='GID')
# 
# # this take time to run, so save it to avoid re-running.
# saveRDS(go_terms, file = here::here("data","intermediate","go_terms.RDS"))
go_terms <- readRDS(file = here::here("data","intermediate","go_terms.RDS"))

# Exclude genes not found in count table --not run--
#go_terms <- go_terms[go_terms$GID %in% rownames(logcpm.norm.counts.combat),]

# gene / go term mapping
gene_go_mapping <- as.data.frame(unique(go_terms %>% select(GID, GO, ONTOLOGY)))
colnames(gene_go_mapping) <- c('gene', 'category', 'ontology')

# go id / term mapping
go_term_id_mapping <- as.data.frame(unique(go_terms[c('GO', 'TERM', 'ONTOLOGY')]))
colnames(go_term_id_mapping) <- c("category", "term", "ontology")

#########################################################################################
# Load KEGG annotations
#########################################################################################

gene_kegg_mapping <- load_kegg_mapping(orgdb, rownames(logcpm.norm.counts.combat),
                                       keytype="GID")

kegg_pathways <- load_kegg_pathways(orgdb, rownames(logcpm.norm.counts.combat),
                                    keytype="GID")


# Rename gene/KEGG mapping columns to be consistent with GO mapping
colnames(gene_kegg_mapping) <- c('gene', 'category')
colnames(kegg_pathways)     <- c('category', 'name', 'class', 'description')

kegg_pathways <- unique(kegg_pathways)


################################################
# GO Enrichment
################################################

# get the number of modules
num_modules <- length(unique(module.colours))

# Create gene lengths vector
gene_lengths <- gene_info$transcript_length
names(gene_lengths) <- gene_info$gene_id

# save the module sizes 
# Data frame of module sizes
module_counts <- c()
for (color in unique(module.colours)) {
  module_counts <- append(module_counts, sum(module.colours == color))
}

# create a mapping from module id to number of genes for later use
module_sizes <- data.frame(module_id=unique(module.colours),
                           num_genes=module_counts)

# Initialize parallelization
cl <- makeCluster(max(1, min(10, detectCores() - 2, na.rm = TRUE)))
registerDoParallel(cl)
message("Performing GO enrichment")

# Check each module for enrichment in GO terms and save result in a list
module_go_enrichment <- foreach(color=unique(module.colours), .packages=c('goseq')) %dopar% {
  set.seed(1)
  # Measure GO enrichment for module
  enriched <- tryCatch({
    # module gene ids
    in_module_geneids <- gene_ids[module.colours == color]
    message(sprintf("[GO enrichment] %s", color))
    
    # T. brucei GO enrichment
    enriched <- test_gene_enrichment(in_module_geneids, gene_ids,
                                     gene_go_mapping, gene_lengths)
    
    # Add descriptions
    enriched <- merge(enriched, go_term_id_mapping, by='category')
  }, error=function(e) {
    # goseq fails in some cases; have not been able to track down cause yet
    # to avoid errors we will just return an empty result set
    warning(sprintf("GO enrichment failed for module %s", color))
    cbind(
      get_enrichment_placeholder(),
      term=numeric(0),
      ontology=numeric(0)
    )
  })
  enriched
}
names(module_go_enrichment) <- unique(module.colours)

# remove any null/empty entries from the results
module_go_enrichment <- module_go_enrichment[!sapply(module_go_enrichment, is.null)]

# unregister cpus
stopCluster(cl)

# save the GO enrichment results
saveRDS(module_go_enrichment, file = here::here("data","intermediate","module_go_enrichment.RDS"))


#------------------------------------
# Print GO enrichment results
#------------------------------------
# temporarily repeat the gene / go term mapping to add 'term' column
gene_go_mapping_tmp <- as.data.frame(unique(go_terms %>% select(GID, GO, TERM, ONTOLOGY)))
colnames(gene_go_mapping_tmp) <- c('gene', 'category', 'term', 'ontology')

gene_info_tmp <- gene_info %>% select(-chromosome, -strand,
                                      - type, -transcript_length)
colnames(gene_info_tmp) <- c("gene","description","transcript_id")

tmp <- cbind(gene_info_tmp, color=module.colours)

#tmp <- cbind(gene=gene_ids, color=module.colours)
gene_mapping <- merge(gene_go_mapping_tmp, tmp, by='gene')
cat(sprintf('- Total enriched modules: %d\n', 
            sum(sapply(module_go_enrichment, nrow) > 0)))

# create tables of the results in this document
print_enrichment_results(module_go_enrichment, module_sizes, 'GO terms',
                         NULL, gene_mapping, 
                         output_dir=here::here("results","tables"),
                         enrichment_type='go',
                         include_gene_lists=FALSE)

enriched_colors_go <- get_enriched_modules(module_go_enrichment)

# Module enrichment status (used in dendrogram plots)
go_enrichment_status   <- as.numeric(module.colours %in% enriched_colors_go)


saveRDS(gene_mapping, file = here::here("data","intermediate","gene_mapping.RDS"))

################################################
# KEGG Enrichment
################################################

# Check each module for enrichment in KEGG terms and save result in a list
cl <- makeCluster(max(1, min(10, detectCores() - 2, na.rm = TRUE)))
registerDoParallel(cl)

message("Performing KEGG enrichment")

# Check each module for enrichment in GO terms and save result in a list
module_kegg_enrichment <- foreach(color=unique(module.colours), .packages=c('goseq')) %dopar% {
  set.seed(1)
  
  # Measure KEGG enrichment for module
  enriched <- tryCatch({
    in_module_geneids <- gene_ids[module.colours == color]
    enriched <- test_gene_enrichment(in_module_geneids, gene_ids,
                                     gene_kegg_mapping, gene_lengths)
    
    enriched <- unique(merge(enriched, kegg_pathways[,c('category','name')],
                             by='category'))
  }, error=function(e) {
    # goseq fails in some cases; have not been able to track down cause yet
    warning(sprintf("KEGG enrichment failed for module %s", color))
    return(get_enrichment_placeholder())
  })
  enriched
}

names(module_kegg_enrichment) <- unique(module.colours)

# remove any null/empty entries from the results
module_kegg_enrichment <- module_kegg_enrichment[!sapply(module_kegg_enrichment, is.null)]

# unregister cpus
stopCluster(cl)

# save the KEGG enrichment results
saveRDS(module_kegg_enrichment, file = here::here("data","intermediate","module_kegg_enrichment.RDS"))


#------------------------------------
# Print KEGG enrichment results
#------------------------------------

cat(sprintf('- Total enriched modules: %d\n', 
            sum(sapply(module_kegg_enrichment, nrow) > 0)))

# create tables of the results in this document --function has bugs for kegg enrichment type--
# print_enrichment_results(module_kegg_enrichment, module_sizes, 
#                          'KEGG pathway',
#                          output_dir= here::here("results","tables"),
#                          enrichment_type='kegg')

enriched_colors_kegg <- get_enriched_modules(module_kegg_enrichment)


# Module enrichment status (used in dendrogram plots)
kegg_enrichment_status <- as.numeric(module.colours %in% enriched_colors_kegg)


# add enrichment information to the Dendrogram
unassigned_modules <- as.numeric(module.colours == 'grey')

png(filename = here::here("results","figures","gene-tree-and-module-enrichment-status.png"), res =1200, 
    type = "cairo", units = 'in', width = 6, height = 6, pointsize = 10)
WGCNA::plotDendroAndColors(gene.tree,
                           cbind(module.colours, go_enrichment_status,
                                 kegg_enrichment_status, unassigned_modules),
                           groupLabels=c(sprintf("Modules\n(n=%s)", num_modules),
                                         #sprintf("Red = upregulated at %s", CONFIG$de_cond2),
                                         "GO enrichment", "KEGG enrichment", "Unassigned"),
                           #cex.colorLabels=cex_color_labels, cex.main=cex_main,
                           # cex.axis=cex_axis, cex.lab=cex_lab,
                           dendroLabels=FALSE,
                           #marAll=c(4,8,6,4),
                           guideHang=0.05)
dev.off()


###############################################################
### Save results
###############################################################

# Save each module's enrichment results as csv files
module_output_dir <- here::here("results","tables","enrichment_results")
if (!dir.exists(module_output_dir)) {
  dir.create(module_output_dir, recursive=TRUE)
}

# GO enrichment 
output_module_enrichment_results(module_go_enrichment, module_output_dir,
                                 'go', go_term_id_mapping)

# KEGG enrichment
output_module_enrichment_results(module_kegg_enrichment, module_output_dir,
                                 'kegg', kegg_pathways %>% select(-description))

#################################################################################
# Save GO enrichment results

# get only the enriched go modules
enriched_module_go_enrichment <- module_go_enrichment[names(module_go_enrichment) %in% enriched_colors_go]

write.xlsx(enriched_module_go_enrichment, 
           file = here::here("results","tables","modules_go_enrichment_results.xlsx"))

####################################################################################
#save KEGG enrichment results

# get only the enriched kegg modules
enriched_module_kegg_enrichment <- module_kegg_enrichment[names(module_kegg_enrichment) %in% enriched_colors_kegg]

write.xlsx(enriched_module_kegg_enrichment, 
           file = here::here("results","tables","modules_kegg_enrichment_results.xlsx"))


#################################################################################
# Create co-expression results data frame
result <- cbind(gene_info, color=module.colours)

# drop unneeded columns
keep_cols <- intersect(c('gene_id', 'color', 'description', 'chromosome',
                         'strand', 'transcript_length'), colnames(result))
result <- tbl_df(result[,keep_cols])

# add expression-related fields
result$expr_variance <- apply(logcpm.norm.counts.combat, 1, var)
result$expr_mean <- apply(logcpm.norm.counts.combat, 1, mean)

# replace colors with numbers
#module.colours <- module_label_mapping$number[match(module.colours, module_label_mapping$color)]

#result$color <- module_label_mapping$number[match(result$color, module_label_mapping$color)]

# Write out an excel/csv file for the result dataframe
openxlsx::write.xlsx(result, file=here::here("results","tables","coexpression_network_genes.xlsx"))


######################################################################################
# save module hub genes domain description
module_hub_genes_description <- gene_info %>% 
                filter(gene_id %in% as.data.frame(module.hub.genes)$module.hub.genes) %>% 
                select(gene_id, description)

module_hub_genes_df <- as.data.frame(module.hub.genes)

module_hub_genes_df <- tibble::rownames_to_column(module_hub_genes_df, var="module")

module_hub_genes_domain_description <- merge(module_hub_genes_df, module_hub_genes_description,
                                             by.x="module.hub.genes", by.y="gene_id")

# reorder and rename columns 
module_hub_genes_domain_description <- module_hub_genes_domain_description[,c(2,1,3)]

module_hub_genes_domain_description <- module_hub_genes_domain_description %>%
  dplyr::rename(hub_gene=module.hub.genes)

openxlsx::write.xlsx(module_hub_genes_domain_description,
           file = here::here("results","tables","hub_genes.xlsx"))

#########################################################################################
# Add hub genes and description to top GO terms in each module

modules_top_over_represented_go_terms <-readxl::read_excel(here::here("results",
                                                                      "tables",
                                                                      "modules_top_over_represented_go_terms.xlsx"))

enriched_modules_hub_genes <- module_hub_genes_domain_description %>% 
              filter(module %in% modules_top_over_represented_go_terms$module)

modules_top_over_represented_go_terms_and_hub_genes <- merge(modules_top_over_represented_go_terms,
                                                             enriched_modules_hub_genes, by="module")

openxlsx::write.xlsx(modules_top_over_represented_go_terms_and_hub_genes,
           file = here::here("results",
                             "tables",
                             "modules_top_over_represented_go_terms_and_hub_genes.xlsx"))
