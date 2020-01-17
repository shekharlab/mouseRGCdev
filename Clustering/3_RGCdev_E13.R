assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

rgcE13 = readRDS("../Objects/rgcE13_190729.rds")
rgcE13 = set.all.ident(rgcE13, id="orig")
DimPlot(rgcE13, reduction.use = "UMAP")
rgcE13 = set.all.ident(rgcE13, id="m")
DimPlot(rgcE13, reduction.use = "UMAP")

# Ligerize
rgcE13@data.info$all = 1
rgcE13 = set.all.ident(rgcE13, id="all")
var.genes = NB.var.genes(rgcE13,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = FALSE, do.text = FALSE)
rgcE13@var.genes = var.genes

library(liger)
rgcE13 = Ligerize(rgcE13, batch_id = "orig", var.genes = rgcE13@var.genes, do.umap = TRUE)
saveRDS(rgcE13, file="../Objects/rgcE13_190729.rds")
rgcE13 = readRDS("../Objects/rgcE13_190729.rds")

rgcE13 = set.all.ident(rgcE13, id="orig")
DimPlot(rgcE13, reduction.use="UMAP")

rgcE13 = set.all.ident(rgcE13, id="m")
DimPlot(rgcE13, reduction.use="UMAP", do.label=TRUE)

p=Perc.pos.by.ident(rgcE13,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25"), max.val.exp = 5, alpha.use=0.7, return.plot = TRUE)
p

# Prune 14,15,16,17,18,19,20
# Cluster 15 was Tmod1, Rhoj, Tesk1, Cnot4
# Cluster 14 was Cdca2, Fgfrl1, Scl9a3r1, Plk3
# Cluster 16 was Atf3, Col4a2, Phgdh, Nid1
rgcE13 = prune.clust(rgcE13, remove.clust = c(14:20) )


ExpMat = Avg.by.ident(rgcE13, features.use = rgcE13@var.genes,return.scaled=FALSE)
a = apply(ExpMat,1, function(x) if(max(x) < 0.3){return(NA)} else {return(mean(sort(x,decreasing=TRUE)[1:7]) / mean(sort(x,decreasing=FALSE)[1:7]))})
a=a[!is.na(a)]
var.genes = names(a)[a>5]
library(lsa)
rgcE13 = buildClusterTree(rgcE13, genes.use=rgcE13@var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "cosine")

# Force merge
rgcE13 = force.merge(rgcE13, merge.clust = c(5,3))

levels(rgcE13@ident) = 1:length(levels(rgcE13@ident))
rgcE13@data.info$m_merge = rgcE13@ident
rgcE13 = find.all.markers(rgcE13, test.use="MAST", max.cells.per.ident = 2000)
saveRDS(rgcE13, file="../Objects/rgcE13_190729.rds")
