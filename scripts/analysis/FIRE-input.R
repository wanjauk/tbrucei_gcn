# write out cluster/module genes and their corresponding module labels for use by 
# FIRE (Finding Informative Regulatory Elements)

fire.clusters <- data.frame(modGenes, module.labels[inModules])
fire.clusters.colours <- data.frame(modGenes, module.labels[inModules], module.colours[inModules])

# sort by module labels; FIRE input should start from 0 in module labels column.
fire.clusters <- fire.clusters[order(fire.clusters$module.labels.inModules.), c(1,2)]
fire.clusters.colours <- fire.clusters.colours[order(fire.clusters.colours$module.labels.inModules.),
                                               c(1,2,3)]

# rename columns
colnames(fire.clusters) <- c("gene", "label")
colnames(fire.clusters.colours) <- c("gene", "label","colour")

write.table(as.data.frame(fire.clusters), file = "../results/tbrucei_FIRE_expression_clusters.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# Also, write out module genes in a text file.
write.table(data.frame(modGenes), 
            file = "../results/tbrucei_module_genes.txt",
            row.names = FALSE, 
            quote = FALSE, 
            col.names = FALSE)

