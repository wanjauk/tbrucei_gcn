#################################################################
# load gene annotations from packages
#################################################################
orgdb <- get("Trypanosoma.brucei.TREU927")

# Fix AnnotationDbi namespace mess
assign('select', dplyr::select, envir=.GlobalEnv)
assign('get',    base::get, envir=.GlobalEnv)

gene_info <- load_parasite_annotations(orgdb, rownames(logcpm.norm.counts.combat),
                                       keytype="GID")

# Get transcript lengths (sum of all exon lengths for each gene)
txdb <- orgdb@txdbSlot
transcript_lengths <- transcriptLengths(txdb)
transcript_lengths <- transcript_lengths[transcript_lengths$gene_id %in%
                                           gene_info$gene_id,]

gene_info[match(transcript_lengths$gene_id, gene_info$gene_id), 
          'transcript_length'] <- transcript_lengths$tx_len
gene_info$transcript_length <- as.numeric(gene_info$transcript_length)

# A gene was excluded in the annotation database as a result of an orphan transcript.
#Add its placeholder.
## ---not run---
gene_info <- rbind(gene_info, data.frame(gene_id='Tb927.4.4663',
                                         chromosome='4',
                                         description=NA, strand=NA, type=NA,
                                         transcript_length=NA))


# Keep only the feature information remaining genes
gene_info <- gene_info[gene_info$gene_id %in% rownames(logcpm.norm.counts.combat),] #WGCNA input counts

# For now, just grab the description for the first transcript
#gene_info <- gene_info[!duplicated(gene_info$gene_id),]

# Gene IDs
gene_ids <- rownames(logcpm.norm.counts.combat)

# gene annotations preview
kable(head(gene_info), caption='Preview of gene annotations.')

#########################################################################################
# load GO terms associated with each parasite gene that were downloaded from Tritrypdb
#########################################################################################
# T. brucei go term annotations file path
#go_term_mapping <- "~/tbrucei_annotation_package/build/3.6.0/TriTrypDB-43_TbruceiTREU927_go_table.txt"
#go_terms <- read.table(go_term_mapping, header = TRUE)

# load go terms from annotation package instead
go_terms <- load_go_terms(orgdb, rownames(logcpm.norm.counts.combat), 
                          keytype='GID')

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


# save required variables obtained earlier in analysis with variable names of the code below
module_colors <- module.colours
gene_tree <-gene.tree
wgcna_input <- logcpm.norm.counts.combat
modules_of_interest <- interesting.modules
num_modules <- length(unique(module_colors))

# save the module sizes 
# Data frame of module sizes
module_counts <- c()
for (color in unique(module_colors)) {
  module_counts <- append(module_counts, sum(module_colors == color))
}

# create a mapping from module id to number of genes for later use
module_sizes <- data.frame(module_id=unique(module_colors),
                           num_genes=module_counts)

# Initialize parallelization
cl <- makeCluster(max(1, min(10, detectCores() - 2, na.rm = TRUE)))
registerDoParallel(cl)
message("Performing GO enrichment")

# Check each module for enrichment in GO terms and save result in a list
module_go_enrichment <- foreach(color=unique(module_colors), .packages=c('goseq')) %dopar% {
  set.seed(1)
  # Measure GO enrichment for module
  enriched <- tryCatch({
    # module gene ids
    in_module_geneids <- gene_ids[module_colors == color]
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
names(module_go_enrichment) <- unique(module_colors)

# remove any null/empty entries from the results
module_go_enrichment <- module_go_enrichment[!sapply(module_go_enrichment, is.null)]

# unregister cpus
stopCluster(cl)

#------------------------------------
# Print enrichment results
#------------------------------------
# temporarily repeat the gene / go term mapping to add 'term' column
gene_go_mapping_tmp <- as.data.frame(unique(go_terms %>% select(GID, GO, TERM, ONTOLOGY)))
colnames(gene_go_mapping_tmp) <- c('gene', 'category', 'term', 'ontology')

gene_info_tmp <- gene_info %>% select(-chromosome, -strand,
                                      - type, -transcript_length)
colnames(gene_info_tmp) <- c("gene","description")

tmp <- cbind(gene_info_tmp, color=module_colors)

#tmp <- cbind(gene=gene_ids, color=module_colors)
gene_mapping <- merge(gene_go_mapping_tmp, tmp, by='gene')
cat(sprintf('- Total enriched modules: %d\n', 
            sum(sapply(module_go_enrichment, nrow) > 0)))

# create tables of the results in this document
print_enrichment_results(module_go_enrichment, module_sizes, 'GO terms',
                         #NULL, gene_mapping, output_dir='../results/printed_enrichment_results', 
                         NULL, gene_mapping, 
                         enrichment_type='go',
                         include_gene_lists=FALSE)

enriched_colors_go <- get_enriched_modules(module_go_enrichment)

# Module enrichment status (used in dendrogram plots)
go_enrichment_status   <- as.numeric(module_colors %in% enriched_colors_go)



################################################
# KEGG Enrichment
################################################


# Check each module for enrichment in KEGG terms and save result in a list
cl <- makeCluster(max(1, min(10, detectCores() - 2, na.rm = TRUE)))
registerDoParallel(cl)

message("Performing KEGG enrichment")

# Check each module for enrichment in GO terms and save result in a list
module_kegg_enrichment <- foreach(color=unique(module_colors), .packages=c('goseq')) %dopar% {
  set.seed(1)
  
  # Measure KEGG enrichment for module
  enriched <- tryCatch({
    in_module_geneids <- gene_ids[module_colors == color]
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

names(module_kegg_enrichment) <- unique(module_colors)

# remove any null/empty entries from the results
module_kegg_enrichment <- module_kegg_enrichment[!sapply(module_kegg_enrichment, is.null)]

# unregister cpus
stopCluster(cl)


cat(sprintf('- Total enriched modules: %d\n', 
            sum(sapply(module_kegg_enrichment, nrow) > 0)))

# create tables of the results in this document
print_enrichment_results(module_kegg_enrichment, module_sizes, 
                         'KEGG pathway',
                         #output_dir='../results/printed_enrichment_results',
                         enrichment_type='kegg')

enriched_colors_kegg <- get_enriched_modules(module_kegg_enrichment)


# Module enrichment status (used in dendrogram plots)
kegg_enrichment_status <- as.numeric(module_colors %in% enriched_colors_kegg)


cat('\n### Dendrogram with annotated modules\n')

unassigned_modules <- as.numeric(module_colors == 'grey')

png(filename = "../figures/figure10_gene-tree-and-module-enrichment-status.png", res =1200, 
    type = "cairo", units = 'in', width = 6, height = 6, pointsize = 10)
WGCNA::plotDendroAndColors(gene_tree,
                           cbind(module_colors, go_enrichment_status,
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


### Module labels mappings
module_label_mapping <- rbind(cbind(modules_of_interest,
                                    seq_along(modules_of_interest)),
                              cbind(remaining_modules, #remaining modules variable?
                                    seq_along(remaining_modules) +
                                      length(modules_of_interest)))

module_label_mapping <- as.data.frame(module_label_mapping)
colnames(module_label_mapping) <- c('color', 'number')


###############################
### Save results
###############################
## Save results

###############################################################################
# Save enrichment results as text files
module_output_dir <- "../results/enrichment_results"
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
# Create result data frame
result <- cbind(gene_info, color=module_colors)

# drop unneeded columns
keep_cols <- intersect(c('gene_id', 'color', 'description', 'chromosome',
                         'strand', 'transcript_length'), colnames(result))
result <- tbl_df(result[,keep_cols])

# add expression-related fields
result$expr_variance <- apply(wgcna_input, 1, var)
result$expr_mean <- apply(wgcna_input, 1, mean)

# replace colors with numbers
#module_colors <- module_label_mapping$number[match(module_colors, module_label_mapping$color)]

#result$color <- module_label_mapping$number[match(result$color, module_label_mapping$color)]

# Write out an excel/csv file for the result dataframe
write.table(result, file="../results/coexpression_network_genes.txt", quote = FALSE, sep = "\t")

write.xlsx(result, file="../results/coexpression_network_genes.xlsx")

#################################################################################
correlation_matrix <- cor(network.counts) # check this further to determine whether relationship is correct

# Module Stats
# main_contrast: variable useful if 'result' variable has Differential expression column specified
# "de_comparisons"
#--this code did not work--
module_stats <- create_module_stats_df(result, correlation_matrix, 
                                       wgcna_input,main_contrast)

# Write module stats to a file
write.csv(module_stats, file="../results/module_stats.csv",row.names=FALSE, quote=FALSE)
#################################################################################
# Save module enrichment results in a dataframe and write out to excel
black <- module_go_enrichment$black
tan <- module_go_enrichment$tan    
brown <- module_go_enrichment$brown
blue <- module_go_enrichment$blue  
turquoise <- module_go_enrichment$turquoise
magenta <- module_go_enrichment$magenta    
darkturquoise <- module_go_enrichment$darkturquoise
green <- module_go_enrichment$green                
red <- module_go_enrichment$red    
pink <- module_go_enrichment$pink
salmon <- module_go_enrichment$salmon
lightyellow <- module_go_enrichment$lightyellow
greenyellow <- module_go_enrichment$greenyellow
purple <- module_go_enrichment$purple

go_enrichment_df_list <- c("black","tan","brown","blue","turquoise","magenta","darkturquoise",
                           "green","red","pink","salmon","lightyellow","greenyellow","purple")

for(name in go_enrichment_df_list){
  write.xlsx(x = get(name), 
             file = "../results/Module_go_enrichment_results.xlsx", 
             sheetName = name, append=TRUE, row.names = FALSE)
}

####################################################################################
#save KEGG enrichment results

red <- module_kegg_enrichment$red    
pink <- module_kegg_enrichment$pink
lightyellow <- module_kegg_enrichment$lightyellow
blue <- module_kegg_enrichment$blue
yellow <- module_kegg_enrichment$yellow
lightcyan <- module_kegg_enrichment$lightcyan
magenta <- module_kegg_enrichment$magenta


kegg_enrichment_df_list <- c("blue","magenta","red","pink","lightyellow","yellow","lightcyan")

for(name in kegg_enrichment_df_list){
  write.xlsx(x = get(name), 
             file = "../results/Module_kegg_enrichment_results.xlsx", 
             sheetName = name, append=TRUE, row.names = FALSE)
}

######################################################################################
# save module hub genes domain description
module_hub_genes_description <- gene_info %>% filter(gene_id %in% as.data.frame(module.hub.genes)$module.hub.genes) %>% select(gene_id, description)

module_hub_genes_df <- as.data.frame(module.hub.genes)

module_hub_genes_df <- tibble::rownames_to_column(module_hub_genes_df, var="module")

module_hub_genes_domain_description <- merge(module_hub_genes_df, module_hub_genes_description,
                                             by.x="module.hub.genes", by.y="gene_id")

# reorder and rename columns 
module_hub_genes_domain_description <- module_hub_genes_domain_description[,c(2,1,3)]

module_hub_genes_domain_description <- module_hub_genes_domain_description %>%
  dplyr::rename(hub_gene=module.hub.genes)

write.xlsx(module_hub_genes_domain_description,
           file = "../results/hub_genes.xlsx",
           #append = TRUE,
           row.names = FALSE)

#########################################################################################
# Add hub genes and description to top GO terms in each module

modules_top_over_represented_go_terms <-readxl::read_excel('../results/modules_top_over_represented_go_terms.xlsx')

enriched_modules_hub_genes <- module_hub_genes_domain_description %>% filter(module %in% modules_top_over_represented_go_terms$module)

modules_top_over_represented_go_terms_and_hub_genes <- merge(modules_top_over_represented_go_terms,
                                                             enriched_modules_hub_genes, by="module")

write.xlsx(modules_top_over_represented_go_terms_and_hub_genes,
           file = "../results/modules_top_over_represented_go_terms_and_hub_genes.xlsx",
           #append = TRUE,
           row.names = FALSE)

######################################################################################
# write to excel DE hub genes expressed through out life cycle stages and their different LogFC in contrasts

contrast_common_de_genes_logFC_hub_genes <- contrast_common_de_genes_logFC %>% filter(gene_id %in% modules_top_over_represented_go_terms_and_hub_genes$hub_genes)

modules_top_over_represented_go_terms_and_hub_genes_logFC <- merge(modules_top_over_represented_go_terms_and_hub_genes,
                                                                   contrast_common_de_genes_logFC_hub_genes,
                                                                   by.x="hub_genes", by.y = "gene_id")

modules_top_over_represented_go_terms_and_hub_genes_logFC <- modules_top_over_represented_go_terms_and_hub_genes_logFC %>% select(-num_total, -num_in_subset)

# reorder columns
modules_top_over_represented_go_terms_and_hub_genes_logFC <- modules_top_over_represented_go_terms_and_hub_genes_logFC[,c(2:6,1,7:9)]

# write output to excel
write.xlsx(modules_top_over_represented_go_terms_and_hub_genes_logFC,
           file = "../results/common_differentially_expressed_hub_genes_logFC_comparison.xlsx",
           #append = TRUE,
           row.names = FALSE)

#######################################################################################
# write hub genes and common DEG across contrasts to output text files for cytoscape analysis

#hub genes
write.table(module_hub_genes_domain_description$hub_genes, file="../results/hub_genes.txt", quote = FALSE, row.names = FALSE)

#common DEGs
write.table(contrast_common_de_genes_logFC$gene_id, file="../results/common_DEG_across_contrasts.txt", quote = FALSE, row.names = FALSE)
