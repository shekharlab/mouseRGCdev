assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Load atlas
load("../../RGConly_MNNCluster_181111.Rdata")

# Ligerize
rgc_atlas@data.info$all = 1
rgc_atlas = set.all.ident(rgc_atlas, id="all")
var.genes = NB.var.genes(rgc_atlas,  set.var.genes = FALSE, num.sd = 1, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0005, rem.mt.rp = FALSE, do.text = FALSE)
rgc_atlas@var.genes = var.genes

rgc_atlas2 = Ligerize(rgc_atlas, batch_id = "batch", var.genes = var.genes, do.umap = TRUE)
rgc_atlas@pca.rot = rgc_atlas2@pca.rot
rgc_atlas@dr = rgc_atlas2@dr
saveRDS(rgc_atlas, file="../Objects/rgc_atlas_190729.rds")
