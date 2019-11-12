assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

earlyRetina = readRDS("../Objects/earlyRetina_E13toP5.rds")
earlyRPC = subsetData(earlyRetina, cells.use = which.cells(earlyRetina, c("P1","P2","P3")))

earlyRPC@data.info$all = 1
earlyRPC = set.all.ident(earlyRPC, id="all")
var.genes = NB.var.genes(earlyRPC,  set.var.genes = FALSE, num.sd = 1, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 30/length(earlyRPC@cell.names), rem.mt.rp = FALSE, do.text = FALSE)
earlyRPC@var.genes = var.genes
earlyRPC = Ligerize(earlyRPC, batch_id = "age", var.genes = earlyRPC@var.genes, do.umap = TRUE)

saveRDS(earlyRPC, file = "../Objects/earlyRPC.rds")

earlyRPC = set.all.ident(earlyRPC, id="m")

DimPlot(earlyRPC, reduction.use = "UMAP", do.label = TRUE)
earlyRPC = prune.clust(earlyRPC, remove.clust = c(7,15,17,25,26))
earlyRPC = set.all.ident(earlyRPC, id="all")
var.genes = NB.var.genes(earlyRPC,  set.var.genes = FALSE, num.sd = 1, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 30/length(earlyRPC@cell.names), rem.mt.rp = FALSE, do.text = FALSE)
earlyRPC@var.genes = var.genes
earlyRPC = Ligerize(earlyRPC, batch_id = "age", var.genes = earlyRPC@var.genes, do.umap = TRUE)
saveRDS(earlyRPC, file = "../Objects/earlyRPC_pruned.rds")

earlyRPC = readRDS("../Objects/earlyRPC_pruned.rds")

# Redo UMAP
library(liger)
earlyRPC =doGraph_clustering(earlyRPC, pcs.use = c(1:15), num.nn=25, method="Louvain", do.jaccard=TRUE, use.reduction = "pca")


earlyRPC = set.all.ident(earlyRPC, id="m")
DimPlot(earlyRPC, reduction.use = "UMAP", do.label = TRUE)
FeaturePlot(earlyRPC,"Neurog2", reduction.use = "UMAP")

library(monocle3)
gene_annotation = data.frame(rownames(earlyRPC@count.data), stringsAsFactors = FALSE)
rownames(gene_annotation) = gene_annotation[,1]
colnames(gene_annotation) = "gene_short_name"


cds <- new_cell_data_set(earlyRPC@count.data,
                         cell_metadata = earlyRPC@data.info,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 100)
reducedDims(cds)$PCA = earlyRPC@dr$pca@cell.embeddings[rownames(reducedDims(cds)$PCA),c(1:15)]
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "m")

earlyRPC =doGraph_clustering(earlyRPC, pcs.use = 1:15, num.nn=25, method="Louvain", do.jaccard=TRUE, use.reduction = "pca")

plot_cells(cds,
           genes=c("H19"),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "m")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "m_prune",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "m_prune",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
