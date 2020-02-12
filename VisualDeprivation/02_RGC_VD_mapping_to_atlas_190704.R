assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

#Read perturbed RGC
load("Objects/rgc_perturb_filt_191122.Rdata")

#Read atlas
rgc_atlas = readRDS("../Developing_RGCs/Objects/rgc_atlas_190911_withP56.rds")


# Find Common Highly variable features
# Training object
rgc_atlas@data.info$all = 1
rgc_atlas = set.all.ident(rgc_atlas, id="all")
var.genes = NB.var.genes(rgc_atlas,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = FALSE, do.text = FALSE)
rgc_atlas = set.all.ident(rgc_atlas, id="mnn_final")

# Testing_object
rgc_perturb@data.info$all = 1
rgc_perturb = set.all.ident(rgc_perturb,id="all")
var.genes2 = NB.var.genes(rgc_perturb,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = FALSE, do.text = FALSE)
rgc_perturb = set.all.ident(rgc_perturb, id="m")

# Variable genes
var.genes_use = intersect(var.genes, var.genes2)


# Use classifier to assign
rgc_atlas = set.all.ident(rgc_atlas, id="mnn_final")
rgc_perturb = set.all.ident(rgc_perturb, id="m")
library(xgboost)
rgc_perturb = XgboostAssign(rgc_atlas, rgc_perturb, var.genes = var.genes_use, do.scale = TRUE, test.label = rgc_perturb@ident, newID = "xgboost0")
rgc_perturb = set.all.ident(rgc_perturb, id = "xgboost0")
cells.ident_new = GraphAssignment(rgc_perturb,knn = 15, reduction.use = "UMAP")
rgc_perturb@data.info$iGraphBoost = cells.ident_new
plotConfusionMatrix(table(rgc_perturb@data.info$m, cells.ident_new), order = "Col")

save(list=c("rgc_perturb"), file="Objects/rgc_perturb_filt_191122.Rdata")
load("Objects/rgc_perturb_filt_191122.Rdata")

# Remove clusters that map very diffusely
DimPlot(rgc_perturb, reduction.use = "UMAP", do.label = TRUE)
rgc_perturb = set.all.ident(rgc_perturb, id="m")
rgc_perturb = find.all.markers(rgc_perturb, thresh.use = 0.3, test.use = "MAST", add.to.slot = TRUE, max.cells.per.ident = 1000)
save(list=c("rgc_perturb"), file="Objects/rgc_perturb_filt_191122.Rdata")

remove.clust = c(26, 28, 29, 33, 35, 36, 38, 40, 41, 42, 43, 44)
A=table(rgc_perturb@data.info$m, rgc_perturb@data.info$batch)
A[remove.clust,]

rgc_perturb = subsetData(rgc_perturb, cells.use = setdiff(rgc_perturb@cell.names, which.cells(rgc_perturb, remove.clust)))
rgc_perturb@data.info$m_new = drop.levels(rgc_perturb@data.info$m)
levels(rgc_perturb@data.info$m_new) = 1:length(levels(rgc_perturb@data.info$m_new))
save(list=c("rgc_perturb"), file="Objects/rgc_perturb_finalfilt_191122.Rdata")

library(xgboost)

rgc_perturb = set.all.ident(rgc_perturb, id="m_new")
pdf("Figs/RGC_perturb_tsne_umap_m_final.pdf", w=9,h=7, useDingbats = FALSE)
DimPlot(rgc_perturb, reduction.use = "UMAP", pt.size=0.3, do.label = TRUE)
DimPlot(rgc_perturb, reduction.use = "tsne", pt.size=0.3, do.label = TRUE)
dev.off()

A = table(rgc_perturb@data.info$iGraphBoostC, rgc_perturb@data.info$m_new)
A = t(scale(t(A), scale=rowSums(A), center=FALSE))
test <- heatmap.2(A, hclustfun = function(x) hclust(x,method="single"))
B=A[test$rowInd, test$colInd]
row.max = apply(B,1,which.max)
row.ord = names(sort(row.max))
pdf("Figs/VD_ConfusionMatrix.pdf", w=8,h=7, useDingbats = FALSE)
plotConfusionMatrix(B[row.ord,], col.low = "white", col.high = "darkblue", ylab.use = "Atlas types", xlab.use = "VD Clusters")
dev.off()


rgc_perturb = set.all.ident(rgc_perturb, id="batch")

pdf("Figs/RGC_perturb_tsne_umap_batch_batch_final.pdf", w=9,h=7, useDingbats=FALSE)
DimPlot(rgc_perturb, reduction.use = "UMAP", pt.size=0.3, do.label = FALSE, cols.use = c("red", "cyan","darkblue"))
DimPlot(rgc_perturb, reduction.use = "tsne", pt.size=0.3, do.label = FALSE, cols.use = c("yellow", "red", "cyan","darkblue"))
dev.off()

rgc_perturb@data.info$iGraphBoostC = factor(paste0("C", as.character(rgc_perturb@data.info$iGraphBoost)), 
                                            levels = paste0("C",c(1:46)))
rgc_perturb = set.all.ident(rgc_perturb, id="iGraphBoostC")
pdf("Figs/RGC_perturb_tsne_umap__final.pdf", w=10,h=7, useDingbats = FALSE)
DimPlot(rgc_perturb, reduction.use = "UMAP", pt.size=0.3, do.label = FALSE)
DimPlot(rgc_perturb, reduction.use = "tsne", pt.size=0.3, do.label = FALSE)
dev.off()

save(list=c("rgc_perturb"), file="Objects/rgc_perturb_finalfilt_191122.Rdata")


