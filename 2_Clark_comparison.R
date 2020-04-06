# Comparison with Clark et al., 2019 data

assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Load Clark Data
library(Matrix)
M = readMM("RawData/Clark_et_al/10x_mouse_retina_development.mtx")
cellIDs = read.csv("RawData/Clark_et_al/10x_mouse_retina_development_phenotype.csv", stringsAsFactors = FALSE)
geneIDs = read.csv("RawData/Clark_et_al/10x_mouse_retina_development_feature.csv", stringsAsFactors = FALSE)


rownames(cellIDs) = as.character(cellIDs$X)
colnames(M) = rownames(cellIDs)
cellIDs = subset(cellIDs, age %in% c("E14","E16","P0"))
M = M[,rownames(cellIDs)]
rownames(M) = geneIDs$gene_short_name
colnames(M) = gsub("P0","P0_r1",colnames(M))
colnames(M) = gsub("E16","E16_r1",colnames(M))

rownames(cellIDs) = gsub("P0","P0_r1",rownames(cellIDs))
rownames(cellIDs) = gsub("E16","E16_r1",rownames(cellIDs))

earlyRetina = readRDS("R_objects/earlyRetina_E13toP5.rds")
earlyRetina = subsetData(earlyRetina, cells.use = which.cells(earlyRetina, c("E14","E16","P0"), id = "age"))

colnames(M) = gsub("_","Clark_", colnames(M))
rownames(cellIDs) = gsub("_","Clark_", rownames(cellIDs))

# Combined count matrix
genes.use = intersect(rownames(earlyRetina@count.data), rownames(M))
M = M[,colSums(M > 0) > 700]
E = cbind(earlyRetina@count.data[genes.use, ], M[genes.use,])

# Create object
earlyRetina_clark=scR(count.data=E,ident.fxn=getStat1)
earlyRetina_clark=setup(earlyRetina_clark,project="earlyRetina_clark",min.cells = 10,min.genes = 700,is.expr=0, threshold.quantile = 0.999)

earlyRetina_clark@data.info$age = "E14"
earlyRetina_clark@data.info[grep("^E16", earlyRetina_clark@cell.names),"age"] = "E16"
earlyRetina_clark@data.info[grep("^P0", earlyRetina_clark@cell.names),"age"] = "P0"
earlyRetina_clark@data.info[grep("^E14Clark", earlyRetina_clark@cell.names),"age"] = "E14Clark"
earlyRetina_clark@data.info[grep("^E16Clark", earlyRetina_clark@cell.names),"age"] = "E16Clark"
earlyRetina_clark@data.info[grep("^P0Clark", earlyRetina_clark@cell.names),"age"] = "P0Clark"

earlyRetina_clark@data.info$age = factor(earlyRetina_clark@data.info$age)

# Broad Annotation
earlyRetina_clark@data.info$broadAnn = "oo"
earlyRetina_clark@data.info[intersect(earlyRetina_clark@cell.names, earlyRetina@cell.names),"broadAnn"] = 
  paste0("ThisStudy_",as.character(earlyRetina@data.info[intersect(earlyRetina_clark@cell.names, earlyRetina@cell.names),"broadAnn"]))

earlyRetina_clark@data.info[intersect(earlyRetina_clark@cell.names, rownames(cellIDs)),"broadAnn"] = 
  paste0("Clark_",as.character(cellIDs[intersect(earlyRetina_clark@cell.names, rownames(cellIDs)),"umap_CellType"]))

earlyRetina_clark@data.info$broadAnn = factor(earlyRetina_clark@data.info$broadAnn)
earlyRetina_clark = set.all.ident(earlyRetina_clark, id="broadAnn")

# Ligerize
earlyRetina_clark@data.info$all = 1
earlyRetina_clark = set.all.ident(earlyRetina_clark, id="all")
var.genes = NB.var.genes(earlyRetina_clark,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0001, rem.mt.rp = FALSE, do.text = FALSE)
earlyRetina_clark@var.genes = var.genes

earlyRetina_clark = Ligerize(earlyRetina_clark, batch_id = "age", var.genes = var.genes, do.umap = TRUE)
saveRDS(earlyRetina_clark, file="R_objects/earlyRetina_clark.rds")

earlyRetina_clark = readRDS("R_objects/earlyRetina_clark.rds")
earlyRetina_clark = set.all.ident(earlyRetina_clark, id="age")



earlyRetina_clark@data.info$sample = "T"
earlyRetina_clark@data.info[grep("Clark", earlyRetina_clark@cell.names),"sample"] = "C"
earlyRetina_clark = set.all.ident(earlyRetina_clark, id="sample")

# Rename 
earlyRetina_clark@data.info$RGCs_only = "Other_ThisStudy"
earlyRetina_clark@data.info[grep("Clark", earlyRetina_clark@cell.names),"RGCs_only"] = "Other_Clark"
earlyRetina_clark@data.info[which.cells(earlyRetina_clark, "ThisStudy_RGC"),"RGCs_only"] = "RGC_ThisStudy"
earlyRetina_clark@data.info[which.cells(earlyRetina_clark, "Clark_Retinal Ganglion Cells"),"RGCs_only"] = "RGC_Clark"
 DimPlot(earlyRetina_clark,group.by = "RGCs_only", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, do.return = TRUE, cols.use = c("gray50","gray70", "blue","cyan"))


p[[1]] = DimPlot(earlyRetina_clark, reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, cols.use = c("black","gray70"), do.return = TRUE)
p[[2]] = FeaturePlot(earlyRetina_clark, "Nefl", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, no.legend = FALSE, do.return = TRUE, cols.use = c("cyan","red"))
p[[2]] = p[[2]][[1]] + theme_classic()
p[[3]] = FeaturePlot(earlyRetina_clark, "Fgf15", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, no.legend = FALSE, do.return = TRUE, cols.use = c("cyan","red"))
p[[3]] = p[[3]][[1]] + theme_classic()
p[[4]] = FeaturePlot(earlyRetina_clark, "Tfap2b", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, no.legend = FALSE, do.return = TRUE, cols.use = c("cyan","red"))
p[[4]] = p[[4]][[1]] + theme_classic()

pdf("Figs/Clark_comparison_UMAP.pdf",w=9,h=7,useDingbats = FALSE)
plot_grid(plotlist = p, nrow=2)
dev.off()

