assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Load atlas
rgc_atlas=readRDS("../R_objects/rgc_atlas_190729.rds")
Count.mat = rgc_atlas@count.data

# load E13
rgcE13 = readRDS("../R_objects/rgcE13_190729.rds")
genes.use = intersect(rownames(Count.mat), rownames(rgcE13@count.data))
Count.mat = cbind(Count.mat[genes.use,], rgcE13@count.data[genes.use,])

# load E14
rgcE14 = readRDS("../R_objects/rgcE14_190729.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcE14@count.data[genes.use,])

# load E16
rgcE16 = readRDS("../R_objects/rgcE16_190729.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcE16@count.data[genes.use,])

# load P0
rgcP0 = readRDS("../R_objects/rgcP0_190729.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcP0@count.data[genes.use,])

# load P5
rgcP5 = readRDS("../R_objects/rgcP5_190729.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcP5@count.data[genes.use,])


# Create object
rgc_dev=scR(count.data=Count.mat,ident.fxn=getStat1)
rgc_dev=setup(rgc_dev,project="earlyRetina",min.cells = 10,min.genes = 700,is.expr=0, threshold.quantile = 0.9995)

rgc_dev@data.info$all = 1
rgc_dev = set.all.ident(rgc_dev, id="all")
var.genes = NB.var.genes(rgc_dev,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 40, do.idents = TRUE, x.low.cutoff = 0.00005, rem.mt.rp = TRUE, do.text = FALSE)
rgc_dev@var.genes = var.genes
var.genes_small = NB.var.genes(rgc_dev,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 0.00005, rem.mt.rp = TRUE, do.text = FALSE)


rgc_dev  = set.all.ident(rgc_dev, id="orig")
rgc_dev@data.info$age = "P56"
rgc_dev@data.info[grep("^E13", rgc_dev@cell.names, value=TRUE), "age"] = "E13"
rgc_dev@data.info[grep("^E14", rgc_dev@cell.names, value=TRUE), "age"] = "E14"
rgc_dev@data.info[grep("^E16", rgc_dev@cell.names, value=TRUE), "age"] = "E16"
rgc_dev@data.info[grep("^P0", rgc_dev@cell.names, value=TRUE), "age"] = "P0"
rgc_dev@data.info[grep("^P5", rgc_dev@cell.names, value=TRUE), "age"] = "P5"
rgc_dev@data.info$age = factor(rgc_dev@data.info$age)


# Batch id
rgc_dev@data.info$batch = as.character(rgc_dev@data.info$orig)
rgc_dev@data.info[grep("^aRGC", rgc_dev@cell.names, value=TRUE), "batch"] = paste0("P56_",as.character(rgc_atlas@data.info[grep("^aRGC", rgc_dev@cell.names, value=TRUE),"batch"]))
rgc_dev@data.info[grep("^P5", rgc_dev@cell.names, value=TRUE), "batch"] = paste0("P5_",as.character(rgcP5@data.info[grep("^P5", rgc_dev@cell.names, value=TRUE),"batch"]))
rgc_dev@data.info$batch = factor(rgc_dev@data.info$batch)

library(liger)
rgc_dev = Ligerize(rgc_dev, batch_id = "age", var.genes = rgc_dev@var.genes, do.umap = TRUE, do.clustering = FALSE)
saveRDS(rgc_dev, file = "../R_objects/rgc_dev_190729.rds")
