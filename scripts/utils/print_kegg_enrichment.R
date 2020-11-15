#' Parses results from goseq enrichment analysis 
#'
#' @author V. Keith Hughitt, \email{khughitt@umd.edu}
#'
#' @param results A dataframe containing one or more results from goseq.
#' @param subset_sizes A dataframe mapping module/DE result names to number 
#'        of genes in that module or DE result.
#' @param annotation_name Type of annotation bring printed
#' @param annotation_mapping An optinal dataframe mapping from annotations
#'        to other useful information such as descriptions.
#' @param gene_mapping An optional dataframe mapping from gene IDs to
#'        annotations.
#' @param output_dir Directory where table output results will be saved to.
#' @param enrichment_type Filename suffix to use to indicate type of enrichment
#'        performed.
#' @param exclude_unclustered Whether or not the uncluster genes (grey) should
#'        be excluded from the results.
#' @param include_gene_lists If true, a list of genes associated with the
#'        enriched annotations will also be printed; requires that
#'        `gene_mapping` also be specified.
#'
#' @return None
print_kegg_enrichment <- function(results, subset_sizes, 
                                     annotation_name='annotations', 
                                     annotation_mapping=NULL,
                                     gene_mapping=NULL,
                                     output_dir=NULL,
                                     enrichment_type="",
                                     exclude_unclustered=TRUE,
                                     include_gene_lists=FALSE,
                                     str_max_width=Inf) {
  # exclude any duplicated pathways to avoid double-counting enriched
  # categories when merging annotation with results
  annotation_mapping <- annotation_mapping[!duplicated(annotation_mapping$category),]
  
  # counters
  total_enriched <- 0
  pvalues <- c()
  
  # show modules with most significant p-values first
  avg_pvals <- sapply(results, function(x) {
    if (nrow(x) > 0) {
      min(median(x$under_represented_pvalue_adj),
          median(x$over_represented_pvalue_adj))
    } else {
      Inf
    }
  })
  results <- results[order(rank(avg_pvals, ties.method='first'))]
  
  # For GO enrichment results, make GO: categories hyperlinks
  if(enrichment_type == 'go') {
    if (!is.null(gene_mapping)) {
      gene_mapping$category <- sprintf("<a href='http://amigo.geneontology.org/amigo/term/%s'>%s</a>",
                                       gene_mapping$category, gene_mapping$category)
    }
    if (!is.null(annotation_mapping)) {
      annotation_mapping$category <- sprintf("<a href='http://amigo.geneontology.org/amigo/term/%s'>%s</a>",
                                             annotation_mapping$category, annotation_mapping$category)
    }
  }
  
  for (result_name in names(results)) {
    result <- results[[result_name]]
    
    # Skip unclustered genes if enabled
    if (exclude_unclustered && result_name == 'grey') {
      next
    }
    
    # Skip entries with no enrichment
    if (nrow(result) == 0) {
      next
    }
    
    # add additional fields
    if (!is.null(annotation_mapping)) {
      result <- merge(result, annotation_mapping, by='category')
    }
    
    # For GO enrichment results, make GO: categories hyperlinks
    if(enrichment_type == 'go') {
      result$category <- sprintf("<a href='http://amigo.geneontology.org/amigo/term/%s'>%s</a>",
                                 result$category, result$category)
    }
    
    # determine number of over- and under-represented categories
    over_rep  <- result[result$over_represented_pvalue_adj < 0.05,]
    under_rep <- result[result$under_represented_pvalue_adj < 0.05,]
    
    # keep track of the number of enriched annotations along with their
    # average adjusted p-value.
    total_enriched <- total_enriched + nrow(over_rep) + nrow(under_rep)
    
    # module size
    num_genes <- subset_sizes$num_genes[subset_sizes$module_id == result_name]
    
    # Print results
    if (nrow(over_rep) > 0 || nrow(under_rep) > 0) {
      # With download links
      if (!is.null(output_dir)) {
        # Table containing enriched annotations for module
        enrichment_filename <- sprintf("%s_%s_enrichment.txt", result_name, enrichment_type)
        enrichment_output <- file.path(output_dir, enrichment_filename)
        
        # Table containing information for genes in the module
        genes_output <- file.path(output_dir, sprintf("%s_genes.txt", result_name))
        
        # Output header with links to tables
        cat(sprintf("\n#### [%s](%s) enrichment (%d [genes](%s))\n", 
                    result_name, enrichment_output, num_genes,
                    genes_output))
      } else {
        # Output header without links to tables
        cat(sprintf("\n#### %s enrichment (%d genes)\n", result_name, num_genes))
      }
    }
    
    # Over-represented terms
    if (nrow(over_rep) > 0) {
      # fields to display
      out <- over_rep %>% select(-over_represented_pvalue,
                                 -under_represented_pvalue,
                                 -under_represented_pvalue_adj)



      cat(sprintf("\n**Over-represented %s:**\n", annotation_name))

      # print
      print(xkable(out %>% dplyr::rename(adj_pval=over_represented_pvalue_adj),
                   str_max_width=str_max_width, row.names=FALSE))
      cat('\n')
      
      # out <- over_rep %>% select(-over_represented_pvalue,
      #                            -under_represented_pvalue,
      #                            -under_represented_pvalue_adj,
      #                            -term.x,
      #                            -ontology.x)
      # 
      # 
      # # Added by Kennedy Mwangi to sort by adjust pvalue and rename columns
      # out <- out[order(out$over_represented_pvalue_adj),]
      # 
      # cat(sprintf("\n**Over-represented %s:**\n", annotation_name))
      # 
      # # print
      # print(xkable(out %>% dplyr::rename(adj_pval=over_represented_pvalue_adj,
      #                                    term=term.y,
      #                                    ontology=ontology.y),
      #              str_max_width=str_max_width, row.names=FALSE))
      # cat('\n')
      
      # Add adjusted pvalues to vector for averaging purposes
      pvalues <- append(pvalues, over_rep$over_represented_pvalue_adj)
      
      # Print specific genes responsible, if requested
      if (include_gene_lists) {
        cat("\n**Genes responsible for enrichment:**\n")             
        gene_list <- gene_mapping %>% filter(category %in% out$category & 
                                               color==result_name)
        gene_list <- gene_list[!duplicated(gene_list),]
        
        print(xkable(gene_list, str_max_width=str_max_width, row.names=FALSE))
        cat('\n')
      }
    }
    
    # Under-represented terms
    if (nrow(under_rep) > 0) {
      # fields to display
      out <- under_rep %>% select(-under_represented_pvalue,
                                  -over_represented_pvalue,
                                  -over_represented_pvalue_adj)

      cat(sprintf("\n**Under-represented %s:**\n", annotation_name))

      print(xkable(out %>% dplyr::rename(adj_pval=under_represented_pvalue_adj),
                   str_max_width=str_max_width, row.names=FALSE))
      cat('\n')


      # out <- under_rep %>% select(-under_represented_pvalue,
      #                             -over_represented_pvalue,
      #                             -over_represented_pvalue_adj,
      #                             -term.x,
      #                             -ontology.x)
      # 
      # # sort output by adjusted pvalue
      # # added by Kennedy Mwangi
      # out <- out[order(out$under_represented_pvalue_adj),]
      # 
      # cat(sprintf("\n**Under-represented %s:**\n", annotation_name))
      # 
      # print(xkable(out %>% dplyr::rename(adj_pval=under_represented_pvalue_adj,
      #                                    term=term.y,
      #                                    ontology=ontology.y),
      #              str_max_width=str_max_width, row.names=FALSE))
      # cat('\n')
      # 
      
      # Add adjusted pvalues to vector for averaging purposes
      pvalues <- append(pvalues, under_rep$under_represented_pvalue_adj)
      
      # Print specific genes responsible, if requested
      if (include_gene_lists) {
        cat("\n**Genes responsible for enrichment:**\n")             
        gene_list <- gene_mapping %>% filter(category %in% out$category & 
                                               color==result_name)
        gene_list <- gene_list[!duplicated(gene_list),]
        
        print(xkable(gene_list, str_max_width=str_max_width, row.names=FALSE))
        cat('\n')
      }
    }
  }
  
  # Summary
  cat(sprintf("\n#### Total enriched %s\n", annotation_name))
  
  if (total_enriched > 0) {
    cat(sprintf("\nTotal: %d (Total -log10(adj.pval) = %f)\n", total_enriched,
                sum(-log10(pmax(pvalues, 1E-10)))))
  } else {
    cat("\nTotal: 0\n")
  }
}