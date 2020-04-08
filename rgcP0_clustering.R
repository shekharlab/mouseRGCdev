assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

rgcP0 = readRDS("../R_objects/rgcP0_190729.rds")
rgcP0 = set.all.ident(rgcP0, id="orig")
DimPlot(rgcP0, reduction.use = "UMAP")
rgcP0 = set.all.ident(rgcP0, id="m")
DimPlot(rgcP0, reduction.use = "UMAP")

# Ligerize
rgcP0@data.info$all = 1
rgcP0 = set.all.ident(rgcP0, id="all")
var.genes = NB.var.genes(rgcP0,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0005, rem.mt.rp = FALSE, do.text = FALSE)
rgcP0@var.genes = var.genes

rgcP0 = Ligerize(rgcP0, batch_id = "orig", var.genes = var.genes, do.umap = TRUE)
saveRDS(rgcP0, file="../R_objects/rgcP0_190729.rds")

rgcP0 = set.all.ident(rgcP0, id="orig")
DimPlot(rgcP0, reduction.use="UMAP")

rgcP0 = set.all.ident(rgcP0, id="m")
DimPlot(rgcP0, reduction.use="UMAP", do.label=TRUE)

p=Perc.pos.by.ident(rgcP0,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25"), max.val.exp = 5, alpha.use=0.7, return.plot = TRUE)
p

a=find.markers(rgcP0, 27, test.use="MAST", max.cells.per.ident = 1000)

# 23 : Clu, Gal, Sst, Cldn11
# 24: Adam7, Jun, Pcdh7, Ramp3, Tln2
# Cluster 27, 
rgcP0 = prune.clust(rgcP0, remove.clust = c(23,27,29))


ExpMat = Avg.by.ident(rgcP0, features.use = rgcP0@var.genes,return.scaled=FALSE)
a = apply(ExpMat,1, function(x) if(max(x) < 0.3){return(NA)} else {return(mean(sort(x,decreasing=TRUE)[1:7]) / mean(sort(x,decreasing=FALSE)[1:7]))})
a=a[!is.na(a)]
var.genes = names(a)[a>5]
library(lsa)
rgcP0 = buildClusterTree(rgcP0, genes.use=rgcP0@var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "cosine")

levels(rgcP0@ident) = 1:length(levels(rgcP0@ident))
rgcP0@data.info$m_merge = rgcP0@ident
rgcP0 = find.all.markers(rgcP0, test.use="MAST", max.cells.per.ident = 2000)
saveRDS(rgcP0, file="../R_objects/rgcP0_190729.rds")

