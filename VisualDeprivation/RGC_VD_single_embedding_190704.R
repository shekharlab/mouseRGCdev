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
load("../RawData/rgc_perturbMouse_201806_10X.Rdata")
Count.mat0 = Count.mat

# Dark rearing
load("../RawData/DarkRearRGC_Sep18_10X.Rdata")
Count.mat0 = cbind(Count.mat0, Count.mat_darkrear)
rm(Count.mat); rm(Count.mat_darkrear)

# bcless
load("../RawData/RGCBCless_countdata_050619.Rdata")
genes.use = intersect(rownames(Count.mat0), rownames(Count.mat_RGCBCless))
Count.mat0 = cbind(Count.mat0[genes.use,], Count.mat_RGCBCless[genes.use,])

# Only use RGCs
load("../RawData/RGC_cellIDs_bcless_rgc_perturb_dr.Rdata")
Count.mat0 = Count.mat0[, rgc_cells]

# Remove "RD1All" samples these are from the "all others"
Count.mat0 = Count.mat0[, grep("RD1All", colnames(Count.mat0), value=TRUE, invert=TRUE)]

# Create object
rgc_perturb=scR(count.data=Count.mat0,ident.fxn=getStat1)
rgc_perturb=setup(rgc_perturb,project="rgc_perturb",min.cells = 10,min.genes = 600,is.expr=0,threshold.quantile = 0.999)


# Variable genes
rgc_perturb@data.info$all = 1
rgc_perturb = set.all.ident(rgc_perturb, id="all")
var.genes_each_no_mt = NB.var.genes(rgc_perturb, set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = TRUE, do.text = TRUE)
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("mt-", var.genes_each_no_mt, value=TRUE))
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("^Rpl", var.genes_each_no_mt, value=TRUE))
var.genes_each_no_mt = setdiff(var.genes_each_no_mt, grep("^Rps", var.genes_each_no_mt, value=TRUE))
rgc_perturb@var.genes = var.genes_each_no_mt

rgc_perturb = set.all.ident(rgc_perturb, id="orig")

# Batch
rgc_perturb@data.info$batch = "DR"
rgc_perturb@data.info[which.cells(rgc_perturb, grep("RD1", levels(rgc_perturb@data.info$orig), value=TRUE)),
                      "batch"] = "RD1"
rgc_perturb@data.info[which.cells(rgc_perturb, grep("RGCBCless", levels(rgc_perturb@data.info$orig), value=TRUE)),
                      "batch"] = "BCless"

rgc_perturb@data.info$batch  = factor(rgc_perturb@data.info$batch,
                                      levels = c("RD1","DR","BCless"))



#PCA + Clustering + Liger
rgc_perturb = Ligerize(rgc_perturb, var.genes = var.genes_each_no_mt, min.cells = 1000, do.umap = TRUE, batch_id = "batch")

dir.create("Figs")

rgc_perturb = set.all.ident(rgc_perturb, id="m")
pdf("Figs/RGC_perturb_tsne_umap.pdf", w=9,h=7, useDingbats = FALSE)
DimPlot(rgc_perturb, reduction.use = "UMAP", pt.size=0.4, do.label = TRUE)
DimPlot(rgc_perturb, reduction.use = "tsne", pt.size=0.4, do.label = TRUE)
dev.off()

rgc_perturb = set.all.ident(rgc_perturb, id="batch")

pdf("Figs/RGC_perturb_tsne_umap_batch.pdf", w=9,h=7, useDingbats=FALSE)
DimPlot(rgc_perturb, reduction.use = "UMAP", pt.size=0.4, do.label = FALSE, cols.use = c("red", "cyan","darkblue"))
DimPlot(rgc_perturb, reduction.use = "tsne", pt.size=0.4, do.label = FALSE, cols.use = c("yellow", "red", "cyan","darkblue"))
dev.off()

dir.create("Objects")
save(list=c("rgc_perturb"), file="Objects/rgc_perturb_190704.Rdata")
rgc_perturb = readRDS("Objects/rgc_perturb_190704.rds")
