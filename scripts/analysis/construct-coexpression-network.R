source(here::here("scripts","analysis","libraries.R"))
source(here::here("scripts","analysis","settings.R"))
logcpm.norm.counts.combat <- readRDS(here::here("data","intermediate","logcpm.norm.counts.combat.RDS"))

########################################################################
## Constructing the network
########################################################################

# obtain the required counts data (WGCNA input)
# WGCNA requires genes to be in columns
network.counts <- t(logcpm.norm.counts.combat)

# determine the soft-thresholding power to use
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(network.counts, powerVector = powers, verbose = 5)

# Plot to determine the soft thresholding power to use.

# Scale-free topology fit index as a function of the soft-thresholding power
png(filename = here::here("results","figures","soft-thresholding-power.png"), res =1200, type = "cairo", units = 'in',
    width = 4, height = 4, pointsize = 10)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");

# The red line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
png(filename = here::here("results","figures","mean-connectivity.png"), res =1200, type = "cairo", units = 'in',
    width = 4, height = 5, pointsize = 10)
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)", 
     ylab="Mean Connectivity", type="n", 
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

# construct adjacency matrix
softpower <- 14
adjacency.matrix <- adjacency(network.counts, power=softpower,
                              type = "signed", corFnc = "cor")

# Turn the adjacency matrix to topologicaal overlap matrix to minimize
# the effects of noise and spurious associations
TOM <- TOMsimilarity(adjacency.matrix, TOMType = "signed")
dissTOM <- 1 - TOM

#set diagonal to NA to remove uninformative correlations
diag(adjacency.matrix) <- NA


# Adjacency matrix heatmap plot / network heatmap of selected genes
heatmap_indices <- sample(nrow(adjacency.matrix), 500) # sub-sample for visualization purposes

png(filename = here::here("results","figures","adjacency-matrix-heatmap.png"), res =1200, type = "cairo", units = 'in',
    width = 5, height = 5, pointsize = 10)
heatmap.2(t(adjacency.matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)
dev.off()



# remove adjacency matrix and TOM to free up memory
rm(adjacency.matrix)
gc()



################################################################
## Detecting co-expression modules in R
################################################################

# view the dendrogram based on hierachical clustering of genes
gene.tree <- flashClust(as.dist(dissTOM), method = "average")

# plot the gene tree
png(filename = here::here("results","figures","gene-tree.png"), res =1200, type = "cairo", units = 'in',
    width = 7, height = 8, pointsize = 10)
#sizeGrWindow(12,9) #open graphical window
plot(gene.tree, xlab="", sub="", main = "Gene clustering based on TOM dissimilarity", 
     labels = FALSE, hang = 0.04)
dev.off()

# identify the modules
module.labels <- cutreeDynamicTree(gene.tree, deepSplit = FALSE, 
                                   minModuleSize = 30)

#view
table(module.labels)

# convert labels to colours
module.colours <- labels2colors(module.labels)

# view
table(module.colours)

# a list of 28 modules

# black          blue         brown          cyan     darkgreen 
# 438           614           547           219           102 
# darkgrey    darkorange       darkred darkturquoise         green 
# 93            80           106           100           528 
# greenyellow          grey        grey60     lightcyan    lightgreen 
# 251            59           191           193           164 
# lightyellow       magenta  midnightblue        orange          pink 
# 129           264           200            86           383 
# purple           red     royalblue        salmon           tan 
# 251           460           127           230           243 
# turquoise         white        yellow 
# 732            61           539 

# visualize the gene tree and TOM matrix together using TOM plot
# if necessary, raise dissTOM to a power to make moderately strong connection more visible in heatmap
diag(dissTOM) <- NA

png(filename = here::here("results","figures","gene-tree-and-dissTOM.png"), res =1200, type = "cairo", units = 'in',
    width = 5, height = 6, pointsize = 10)
TOMplot(dissTOM, gene.tree, as.character(module.colours))
dev.off()

# remove matrix to free memory
rm(dissTOM)
gc()


# plot gene dendrogram
png(filename = here::here("results","figures","gene-tree-and-colours.png"), res =1200, type = "cairo", units = 'in',
    width = 6, height = 6, pointsize = 10)
#sizeGrWindow(8,6) #open graphical window
plotDendroAndColors(gene.tree, module.colours, "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colours")
dev.off()

# get hub genes
# choose power 4: https://support.bioconductor.org/p/46342/
module.hub.genes <- chooseTopHubInEachModule(network.counts, module.colours, 
                                             power = 4,type = "signed")

# # A list of module hub genes
module.hub.genes
# black             blue            brown             cyan 
# "Tb927.7.1790"   "Tb927.8.3620"  "Tb927.11.1570"   "Tb927.2.5530" 
# darkgreen         darkgrey       darkorange          darkred 
# "Tb927.11.14020"   "Tb927.9.9450"    "Tb927.8.710"   "Tb927.2.5270" 
# darkturquoise            green      greenyellow           grey60 
# "Tb927.8.6650" "Tb927.10.13790"    "Tb927.7.920"  "Tb927.10.3820" 
# lightcyan       lightgreen      lightyellow          magenta 
# "Tb927.11.6440"   "Tb927.10.720"  "Tb927.10.2560"   "Tb927.9.6290" 
# midnightblue           orange             pink           purple 
# "Tb927.11.6640"   "Tb927.9.2520"  "Tb927.10.6200"    "Tb927.1.600" 
# red        royalblue           salmon              tan 
# "Tb927.7.6920" "Tb927.10.15680"  "Tb927.11.1450"   "Tb927.3.2930" 
# turquoise            white           yellow 
# "Tb927.9.15630"   "Tb927.8.7980"   "Tb927.1.3550"


##############################################################################
## Network export to cytoscape
##############################################################################

# select modules of interest
all.modules <- c('black', 'cyan', 'grey','brown',
                         'midnightblue','blue','darkgreen','tan', 'darkgrey','darkorange',
                         'darkred','darkturquoise','green','greenyellow','grey60','lightcyan',
                         'lightgreen','lightyellow','magenta','orange','pink','purple','red',
                         'royalblue','salmon','turquoise','white','yellow') # all module colours

# enriched modules
enriched.modules <- c("black","tan","brown","blue","turquoise","magenta","darkturquoise",
                      "green","red","pink","salmon","lightyellow","purple","greenyellow")

# obtain gene ids
gene.ids <- rownames(logcpm.norm.counts.combat)

# select module genes
inModules <- is.finite(match(module.colours, all.modules)) # whole network modules
#inModules <- is.finite(match(module.colours, enriched.modules)) # enriched modules

modGenes <- gene.ids[inModules]

# select the corresponding dissTOM based on module genes
modTOM <- TOM[inModules, inModules]
dimnames(modTOM) <- list(modGenes, modGenes)

# Export the network into edge and node list files Cytoscape can read
exportNetworkToCytoscape(modTOM,
                         edgeFile = here::here("data","intermediate","cytoscape_input_files",
                                               "CytoscapeInput-edges_whole_network_thresh0.3.txt"),
                         nodeFile = here::here("data","intermediate","cytoscape_input_files",
                                               "CytoscapeInput-nodes_whole_network_thresh0.3.txt"),

                         weighted = TRUE,
                         threshold = 0.3,
                         nodeNames = modGenes,
                         nodeAttr = module.colours[inModules]);

# network modules
# create a dataframe with node attributes
enriched.module.colours <- module.colours[inModules] #get enriched module colours from module.colours
node.attributes <- cbind(modGenes, module=module.colours) # get node attr. for whole network
# node.attributes <- cbind(modGenes, module=enriched.module.colours) # node atrr. for enriched modules

node.attributes <- as.data.frame(node.attributes)

# Add RGB versions of colour modules
node.attributes$colourRGB <- col2hex(node.attributes$module)

# write out a node attributes files with hexadecimal colour names for module genes
write.table(node.attributes, 
            file = here::here("data","intermediate","cytoscape_input_files",
                              "Cytoscape_node_attributes_whole_network_thresh0.3.txt"),
            row.names = FALSE, 
            quote = FALSE, sep = "\t")

########################################
## FIRE Motif Prediction Input
########################################
# write out cluster/module genes and their corresponding module labels for use by 
# FIRE (Finding Informative Regulatory Elements)

fire.clusters.colours <- data.frame(modGenes, module.labels[inModules], module.colours[inModules])

# sort by module labels; FIRE input should start from 0 in module labels column.
fire.clusters.colours <- fire.clusters.colours[order(fire.clusters.colours$module.labels.inModules.),
                                               c(1,2,3)]

# rename columns
colnames(fire.clusters.colours) <- c("gene", "label","colour")

write.table(as.data.frame(fire.clusters.colours), file = here::here("data","intermediate","tbrucei_FIRE_expression_clusters_and_colour.txt"), 
            quote = FALSE, row.names = FALSE, sep = "\t")

# Also, write out module genes in a text file.
write.table(data.frame(modGenes), 
            file = here::here("data","intermediate","tbrucei_module_genes.txt"),
            row.names = FALSE, 
            quote = FALSE, 
            col.names = FALSE)


saveRDS(gene.tree, file = here::here("data","intermediate","gene.tree.RDS"))
saveRDS(all.modules, file = here::here("data","intermediate","all.modules.RDS"))
saveRDS(module.colours, file = here::here("data","intermediate","module.colours.RDS"))
saveRDS(module.hub.genes, file = here::here("data","intermediate","module.hub.genes.RDS"))
