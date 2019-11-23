assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Load datasets
# rgc_perturb
load("../RawData/RD1Mouse_201806_10X.Rdata")
Count.mat0 = Count.mat

# Dark rearing
load("../RawData/DarkRearRGC_Sep18_10X.Rdata")
Count.mat0 = cbind(Count.mat0, Count.mat_darkrear)
rm(Count.mat); rm(Count.mat_darkrear)

# bcless
load("../RawData/RGCBCless_countdata_050619.Rdata")
genes.use = intersect(rownames(Count.mat0), rownames(Count.mat_RGCBCless))
Count.mat0 = cbind(Count.mat0[genes.use,], Count.mat_RGCBCless[genes.use,])


# Remove "RD1All" samples these are from the "all others"
Count.mat0 = Count.mat0[, grep("RD1All", colnames(Count.mat0), value=TRUE, invert=TRUE)]

# Create object
rgc_perturb=scR(count.data=Count.mat0,ident.fxn=getStat1)
rgc_perturb=setup(rgc_perturb,project="rgc_perturb",min.cells = 10,min.genes = 600,is.expr=0,threshold.quantile = 0.999)

# Sample Origin
rgc_perturb@data.info$batch = "DR"
rgc_perturb@data.info[which.cells(rgc_perturb, grep("RD1", levels(rgc_perturb@data.info$orig), value=TRUE)),
                      "batch"] = "RD1"
rgc_perturb@data.info[which.cells(rgc_perturb, grep("RGCBCless", levels(rgc_perturb@data.info$orig), value=TRUE)), "batch"] = "BCless"

rgc_perturb@data.info$batch  = factor(rgc_perturb@data.info$batch,
                                      levels = c("RD1","DR","BCless"))


########
### First round of clustering
#########
# Variable genes
rgc_perturb@data.info$all = 1
rgc_perturb = set.all.ident(rgc_perturb, id="all")
var.genes_each_no_mt = NB.var.genes(rgc_perturb, set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = TRUE, do.text = TRUE)
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("mt-", var.genes_each_no_mt, value=TRUE))
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("^Rpl", var.genes_each_no_mt, value=TRUE))
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("^Rps", var.genes_each_no_mt, value=TRUE))
rgc_perturb@var.genes = var.genes_each_no_mt
rgc_perturb = set.all.ident(rgc_perturb, id="orig")

#PCA + Clustering + Liger
rgc_perturb = Ligerize(rgc_perturb, var.genes = var.genes_each_no_mt, min.cells = 1000, do.umap = TRUE, batch_id = "batch")

dir.create("Objects")
save(list=c("rgc_perturb"), file="Objects/rgc_perturb_unfilt_191122.Rdata")

ExpMat = Avg.by.ident(rgc_perturb, features.use = rgc_perturb@var.genes,return.scaled=FALSE)
a = apply(ExpMat,1, function(x) if(max(x) < 0.3){return(NA)} else {return(mean(sort(x,decreasing=TRUE)[1:7]) / mean(sort(x,decreasing=FALSE)[1:7]))})
a=a[!is.na(a)]
var.genes = names(a)[a>5]
library(lsa)
rgc_perturb = buildClusterTree(rgc_perturb, genes.use=var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "cosine")

plot_dendro_withdotplot(rgc_perturb,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25","Prkca","Vsx2","Cabp5","Scgn","Rho","Nrl","Pdc"), max.val.exp = 5, alpha.use=0.7)

# Make Broad annotations
rgc_perturb@data.info$cell_class = "RGC"
rgc_perturb@data.info[which.cells(rgc_perturb, c(40,29,2,34, 28, 12, 38)),"cell_class"] = "AC"
rgc_perturb@data.info[which.cells(rgc_perturb, c(8,41)),"cell_class"] = "PR"
rgc_perturb@data.info[which.cells(rgc_perturb, c(36,39)),"cell_class"] = "Other"

library(Hmisc)
pdf("Figs/ClassDistribution_VD_prefilt.pdf",w=5,h=5,useDingbats = FALSE)
cluster.dist.by.ident(rgc_perturb, ident1.use="cell_class", ident2.use = "batch")
dev.off()
save(list=c("rgc_perturb"), file="Objects/rgc_perturb_unfilt_191122.Rdata")


# Extract RGCs
rgc_perturb = set.all.ident(rgc_perturb, id="cell_class")
rgc_perturb = subsetData(rgc_perturb, cells.use = which.cells(rgc_perturb,"RGC"))




########
### Second round of clustering
#########
# Variable genes
rgc_perturb@data.info$all = 1
rgc_perturb = set.all.ident(rgc_perturb, id="all")
var.genes_each_no_mt = NB.var.genes(rgc_perturb, set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = TRUE, do.text = TRUE)
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("mt-", var.genes_each_no_mt, value=TRUE))
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("^Rpl", var.genes_each_no_mt, value=TRUE))
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("^Rps", var.genes_each_no_mt, value=TRUE))
rgc_perturb@var.genes = var.genes_each_no_mt
rgc_perturb = set.all.ident(rgc_perturb, id="orig")

#PCA + Clustering + Liger
rgc_perturb = Ligerize(rgc_perturb, var.genes = var.genes_each_no_mt, min.cells = 1000, do.umap = TRUE, batch_id = "batch")

dir.create("Objects")
save(list=c("rgc_perturb"), file="Objects/rgc_perturb_filt_191122.Rdata")

ExpMat = Avg.by.ident(rgc_perturb, features.use = rgc_perturb@var.genes,return.scaled=FALSE)
a = apply(ExpMat,1, function(x) if(max(x) < 0.3){return(NA)} else {return(mean(sort(x,decreasing=TRUE)[1:7]) / mean(sort(x,decreasing=FALSE)[1:7]))})
a=a[!is.na(a)]
var.genes = names(a)[a>5]
library(lsa)
rgc_perturb = buildClusterTree(rgc_perturb, genes.use=var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "cosine")

plot_dendro_withdotplot(rgc_perturb,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25","Prkca","Vsx2","Cabp5","Scgn","Rho","Nrl","Pdc"), max.val.exp = 5, alpha.use=0.7)

save(list=c("rgc_perturb"), file="Objects/rgc_perturb_filt_191122.Rdata")

dir.create("Figs")

rgc_perturb = set.all.ident(rgc_perturb, id="m")
pdf("Figs/RGC_perturb_tsne_umap2.pdf", w=9,h=7, useDingbats = FALSE)
DimPlot(rgc_perturb, reduction.use = "UMAP", pt.size=0.4, do.label = TRUE)
DimPlot(rgc_perturb, reduction.use = "tsne", pt.size=0.4, do.label = TRUE)
dev.off()

rgc_perturb = set.all.ident(rgc_perturb, id="batch")

pdf("Figs/RGC_perturb_tsne_umap_batch2.pdf", w=9,h=7, useDingbats=FALSE)
DimPlot(rgc_perturb, reduction.use = "UMAP", pt.size=0.4, do.label = FALSE, cols.use = c("red", "cyan","darkblue"))
DimPlot(rgc_perturb, reduction.use = "tsne", pt.size=0.4, do.label = FALSE, cols.use = c("yellow", "red", "cyan","darkblue"))
dev.off()
