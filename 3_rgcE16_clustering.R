assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

rgcE16 = readRDS("../R_objects/rgcE16_190729.rds")
rgcE16 = set.all.ident(rgcE16, id="orig")
DimPlot(rgcE16, reduction.use = "UMAP")
rgcE16 = set.all.ident(rgcE16, id="m")
DimPlot(rgcE16, reduction.use = "UMAP")

# Ligerize
rgcE16@data.info$all = 1
rgcE16 = set.all.ident(rgcE16, id="all")
var.genes = NB.var.genes(rgcE16,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0005, rem.mt.rp = FALSE, do.text = FALSE)
rgcE16@var.genes = var.genes

rgcE16 = Ligerize(rgcE16, batch_id = "orig", var.genes = var.genes, do.umap = TRUE)
saveRDS(rgcE16, file="../R_objects/rgcE16_190729.rds")

rgcE16 = set.all.ident(rgcE16, id="orig")
DimPlot(rgcE16, reduction.use="UMAP")

rgcE16 = set.all.ident(rgcE16, id="m")
DimPlot(rgcE16, reduction.use="UMAP", do.label=TRUE)

p=Perc.pos.by.ident(rgcE16,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25"), max.val.exp = 5, alpha.use=0.7, return.plot = TRUE)
p

a=find.markers(rgcE16, 16, test.use="MAST", max.cells.per.ident = 1000)
rgcE16 = prune.clust(rgcE16, remove.clust = c(11,18,19,22,23,25))


ExpMat = Avg.by.ident(rgcE16, features.use = rgcE16@var.genes,return.scaled=FALSE)
a = apply(ExpMat,1, function(x) if(max(x) < 0.3){return(NA)} else {return(mean(sort(x,decreasing=TRUE)[1:7]) / mean(sort(x,decreasing=FALSE)[1:7]))})
a=a[!is.na(a)]
var.genes = names(a)[a>5]
library(lsa)
rgcE16 = buildClusterTree(rgcE16, genes.use=rgcE16@var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "cosine")

levels(rgcE16@ident) = 1:length(levels(rgcE16@ident))
rgcE16@data.info$m_merge = rgcE16@ident
rgcE16 = find.all.markers(rgcE16, test.use="MAST", max.cells.per.ident = 2000)
saveRDS(rgcE16, file="../R_objects/rgcE16_190729.rds")

