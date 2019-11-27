assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

#### we will use the rgc_dev object as it is normalized appropriately
rgc_dev = readRDS("../Developing_RGCs/Objects/rgc_dev_190910.rds")

#### Atlas
rgc_atlas = readRDS("../Developing_RGCs/Objects/rgc_atlas_190911_withP56.rds")
orig_atlas_expmat = rgc_atlas@data
rgc_atlas@data = rgc_dev@data[,rgc_atlas@cell.names]

#### P5
rgcP5 = readRDS("../Developing_RGCs/Objects/rgcP5_190911_withP56.rds")
orig_P5_expmat = rgcP5@data
rgcP5@data = rgc_dev@data[,rgcP5@cell.names]

#### Assign P5 RGCs an atlas identity using iGraphBoost
# Find Common Highly variable features
# Training object
rgc_atlas@data.info$all = 1
rgc_atlas = set.all.ident(rgc_atlas, id="all")
var.genes = NB.var.genes(rgc_atlas,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = FALSE, do.text = FALSE)
rgc_atlas = set.all.ident(rgc_atlas, id="mnn_final")

# Testing_object
rgcP5@data.info$all = 1
rgcP5 = set.all.ident(rgcP5,id="all")
var.genes2 = NB.var.genes(rgcP5,  set.var.genes = FALSE, num.sd = 0.7, x.high.cutoff = 10, do.idents = TRUE, x.low.cutoff = 0.001, rem.mt.rp = FALSE, do.text = FALSE)
rgcP5 = set.all.ident(rgcP5, id="mnn_final")

# Variable genes
var.genes_use = intersect(var.genes, var.genes2)


# Use classifier to assign
rgc_atlas = set.all.ident(rgc_atlas, id="mnn_final")
rgcP5 = set.all.ident(rgcP5, id="mnn_final")
source("iGraphBoost_scripts.R")
library(xgboost)
rgcP5 = XgboostAssign(rgc_atlas, rgcP5, var.genes = var.genes_use, do.scale = TRUE, test.label = rgcP5@ident, newID = "xgboost0")
rgcP5 = set.all.ident(rgcP5, id = "xgboost0")
cells.ident_new = GraphAssignment(rgcP5,knn = 15, reduction.use = "UMAP")
rgcP5@data.info$iGraphBoost = cells.ident_new
plotConfusionMatrix(table(rgcP5@data.info$mnn_final, cells.ident_new), order = "Col")

# OT assignments
OT = rgcP5@data.info[,grep("fate$", colnames(rgcP5@data.info))]
OT = OT[,paste0(levels(rgc_atlas@data.info$new_types),"_fate")]
rgcP5@data.info$OT_assignment = apply(OT,1,function(x) if (max(x) > 0.1){which.max(x)} else {46})
plotConfusionMatrix(table(rgcP5@data.info$mnn_final, rgcP5@data.info$OT_assignment), order="Col")
plotConfusionMatrix(table(rgcP5@data.info$iGraphBoost, rgcP5@data.info$OT_assignment), order="Col")
rgcP5@data.info$OT_assignment = factor(rgcP5@data.info$OT_assignment)

# Resave the rgcP5 object
rgcP5@data = orig_P5_expmat
saveRDS(rgcP5, file = "../Developing_RGCs/Objects/rgcP5_190911_withP56.rds")
