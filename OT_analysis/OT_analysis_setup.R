assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Load atlas
rgc_atlas=readRDS("../Developing_RGCs/Objects/rgc_atlas_190830_withP56.rds")
Count.mat = rgc_atlas@count.data

# load E13
rgcE13 = readRDS("../Developing_RGCs/Objects/rgcE13_190830_withP56.rds")
genes.use = intersect(rownames(Count.mat), rownames(rgcE13@count.data))
Count.mat = cbind(Count.mat[genes.use,], rgcE13@count.data[genes.use,])

# load E14
rgcE14 = readRDS("../Developing_RGCs/Objects/rgcE14_190830_withP56.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcE14@count.data[genes.use,])

# load E16
rgcE16 = readRDS("../Developing_RGCs/Objects/rgcE16_190830_withP56.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcE16@count.data[genes.use,])

# load P0
rgcP0 = readRDS("../Developing_RGCs/Objects/rgcP0_190830_withP56.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcP0@count.data[genes.use,])

# load P5
rgcP5 = readRDS("../Developing_RGCs/Objects/rgcP5_190830_withP56.rds")
Count.mat = cbind(Count.mat[genes.use,], rgcP5@count.data[genes.use,])


# Create object
rgc_dev=scR(count.data=Count.mat,ident.fxn=getStat1)
rgc_dev=setup(rgc_dev,project="earlyRetina",min.cells = 10,min.genes = 700,is.expr=0, threshold.quantile = 0.9995)

rgc_dev@data.info$all = 1
rgc_dev = set.all.ident(rgc_dev, id="all")
var.genes = NB.var.genes(rgc_dev,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 40, do.idents = TRUE, x.low.cutoff = 0.00005, rem.mt.rp = TRUE, do.text = FALSE)
rgc_dev@var.genes = var.genes
var.genes_small = NB.var.genes(rgc_dev,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 0.00005, rem.mt.rp = TRUE, do.text = FALSE)
saveRDS(rgc_dev, file = "../Developing_RGCs/Objects/rgc_dev_190910.rds")
rgc_dev = readRDS("../Developing_RGCs/Objects/rgc_dev_190910.rds")

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
saveRDS(rgc_dev, file = "../Developing_RGCs/Objects/rgc_dev_190910.rds")


rgc_dev = readRDS("../Developing_RGCs/Objects/rgc_dev_190816.rds")

dir.create("OT_uncorr_vargenes_eps005")

# Var genes small
writeMM(rgc_dev@data[var.genes_small,], "OT_uncorr_vargenes_eps005/RGC_batch_uncorr_190910_small.mtx")
write(colnames(rgc_dev@data), file = "OT_uncorr_vargenes_eps005/RGC_batch_uncorr_190910_small.barcodes.txt", ncolumns=1)
write(var.genes_small, file = "OT_uncorr_vargenes_eps005/RGC_batch_uncorr_190910_small.genes.txt", ncolumns=1)


cell.ids = colnames(rgc_dev@data)
df = data.frame(id=cell.ids, day=0)
rownames(df) = cell.ids
df[grep("^E14",cell.ids,value=TRUE),"day"] = 1
df[grep("^E16",cell.ids,value=TRUE),"day"] = 3
df[grep("^P0",cell.ids,value=TRUE),"day"] = 6
df[grep("^P5",cell.ids,value=TRUE),"day"] = 11
df[grep("^aRGC",cell.ids,value=TRUE),"day"] = 16

write.table(df, file = "OT_uncorr_vargenes_eps005/cell_day.txt", row.names = FALSE, sep="\t", quote=FALSE)

# Write cell sets
objects = c("rgc_atlas","rgcP5", "rgcP0", "rgcE16", "rgcE14","rgcE13")
ages = c("P56","P5","P0","E16","E14","E13")
sink("OT_uncorr_vargenes_eps005/cell_sets.gmt")
#library(xgboost)
for (n in c(1:6)){
  eval(parse(text=paste0("object = ", objects[n])))
  for (i in levels(object@ident)){
    cells.use = which.cells(object, i)
    if (n != 1){
      k = paste0("C",i)
    } else {
      k = i
    }
    set_name = paste0(ages[n], "_", k)
    cells_set = paste0(cells.use, collapse = "\t")
    
    str_to_write = paste0(set_name,"\t NA \t", cells_set)
    cat(str_to_write)
    cat("\n")
  }
}
sink()
