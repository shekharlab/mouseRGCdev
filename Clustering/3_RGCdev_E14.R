assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

rgcE14 = readRDS("../Objects/rgcE14_190729.rds")
rgcE14 = set.all.ident(rgcE14, id="orig")
DimPlot(rgcE14, reduction.use = "UMAP")
rgcE14 = set.all.ident(rgcE14, id="m")
DimPlot(rgcE14, reduction.use = "UMAP")

# Ligerize
rgcE14@data.info$all = 1
rgcE14 = set.all.ident(rgcE14, id="all")
var.genes = NB.var.genes(rgcE14,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0001, rem.mt.rp = FALSE, do.text = FALSE)
rgcE14@var.genes = var.genes

rgcE14 = Ligerize(rgcE14, batch_id = "orig", var.genes = var.genes, do.umap = TRUE)
saveRDS(rgcE14, file="../Objects/rgcE14_190729.rds")
rgcE14 = readRDS("../Objects/rgcE14_190729.rds")

rgcE14 = set.all.ident(rgcE14, id="orig")
DimPlot(rgcE14, reduction.use="UMAP")

rgcE14 = set.all.ident(rgcE14, id="m")
DimPlot(rgcE14, reduction.use="UMAP", do.label=TRUE)

plot_dendro_withdotplot(rgcE14,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25"), max.val.exp = 5, alpha.use=0.7)

# Prune 13, 21,22, 24, 26,20
rgcE14 = prune.clust(rgcE14, remove.clust = c(13,20,21,22,24,26,27) )


ExpMat = Avg.by.ident(rgcE14, features.use = rgcE14@var.genes,return.scaled=FALSE)
a = apply(ExpMat,1, function(x) if(max(x) < 0.3){return(NA)} else {return(mean(sort(x,decreasing=TRUE)[1:7]) / mean(sort(x,decreasing=FALSE)[1:7]))})
a=a[!is.na(a)]
var.genes = names(a)[a>5]
library(lsa)
rgcE14 = buildClusterTree(rgcE14, genes.use=rgcE14@var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "cosine")

# Force merge
rgcE14 = force.merge(rgcE14, merge.clust = c(1,2,4,11,16,20))
levels(rgcE14@ident) = 1:length(levels(rgcE14@ident))
rgcE14@data.info$m_merge = rgcE14@ident
rgcE14 = find.all.markers(rgcE14, test.use="MAST", max.cells.per.ident = 2000)
saveRDS(rgcE14, file="../Objects/rgcE14_190729.rds")

