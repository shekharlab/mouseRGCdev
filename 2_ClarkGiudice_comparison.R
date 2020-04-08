# Comparison with Clark et al., 2019 and Lo Giudice, 2019

assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")
library(liger)

# Load Lio Giudice 
Counts_LG = read.csv("RawData/GSE122466_Merged5347cells_RAW.csv", stringsAsFactors = FALSE, sep=";")
rownames(Counts_LG) = Counts_LG[,1]
Counts_LG = Counts_LG[,-1]
Counts_LG = Matrix(as.matrix(Counts_LG), sparse=TRUE)

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
earlyRetina = set.all.ident(earlyRetina, id="age")
cells.use = subsample_by_ident(earlyRetina@ident, max.cells.per.ident = 12000, downsample.frac = 1)
earlyRetina = subsetData(earlyRetina, cells.use=cells.use)

# Add tags to colnames
colnames(M) = gsub("_","Clark_", colnames(M))
rownames(cellIDs) = gsub("_","Clark_", rownames(cellIDs))
colnames(Counts_LG) = paste0("E15.5LoGiudice_", colnames(Counts_LG))

# Combined count matrix
genes.use = intersect(intersect(rownames(earlyRetina@count.data), rownames(M)), rownames(Counts_LG))
M = M[,colSums(M > 0) > 700]
E = cbind(earlyRetina@count.data[genes.use, cells.use ], M[genes.use,])
Counts_LG = Counts_LG[,colSums(Counts_LG > 0) > 700]
E = cbind(E, Counts_LG[genes.use,])


# Create object
earlyRetina_all=scR(count.data=E,ident.fxn=getStat1)
earlyRetina_all=setup(earlyRetina_all,project="earlyRetina_all",min.cells = 10,min.genes = 700,is.expr=0, threshold.quantile = 0.999)

earlyRetina_all@data.info$age = "E14"
earlyRetina_all@data.info[grep("^E16", earlyRetina_all@cell.names),"age"] = "E16"
earlyRetina_all@data.info[grep("^E15.5Lo", earlyRetina_all@cell.names),"age"] = "E15.5LoGiudice"
earlyRetina_all@data.info[grep("^P0", earlyRetina_all@cell.names),"age"] = "P0"
earlyRetina_all@data.info[grep("^E14Clark", earlyRetina_all@cell.names),"age"] = "E14Clark"
earlyRetina_all@data.info[grep("^E16Clark", earlyRetina_all@cell.names),"age"] = "E16Clark"
earlyRetina_all@data.info[grep("^P0Clark", earlyRetina_all@cell.names),"age"] = "P0Clark"

earlyRetina_all@data.info$age = factor(earlyRetina_all@data.info$age)
earlyRetina_all = set.all.ident(earlyRetina_all, id="age")
VlnPlot(earlyRetina_all, "nGene", x.lab.rot=TRUE)

# Broad Annotation
earlyRetina_all@data.info$broadAnn = "oo"
earlyRetina_all@data.info[intersect(earlyRetina_all@cell.names, earlyRetina@cell.names),"broadAnn"] = 
  paste0("ThisStudy_",as.character(earlyRetina@data.info[intersect(earlyRetina_all@cell.names, earlyRetina@cell.names),"broadAnn"]))

earlyRetina_all@data.info[intersect(earlyRetina_all@cell.names, rownames(cellIDs)),"broadAnn"] = 
  paste0("Clark_",as.character(cellIDs[intersect(earlyRetina_all@cell.names, rownames(cellIDs)),"umap_CellType"]))

earlyRetina_all@data.info[intersect(earlyRetina_all@cell.names, colnames(Counts_LG)),"broadAnn"] = "LoGiudice"

earlyRetina_all@data.info$broadAnn = factor(earlyRetina_all@data.info$broadAnn)
earlyRetina_all = set.all.ident(earlyRetina_all, id="broadAnn")

# Ligerize
earlyRetina_all@data.info$all = 1
earlyRetina_all = set.all.ident(earlyRetina_all, id="all")
var.genes = NB.var.genes(earlyRetina_all,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0001, rem.mt.rp = FALSE, do.text = FALSE)
earlyRetina_all@var.genes = var.genes

earlyRetina_all = Ligerize(earlyRetina_all, batch_id = "age", var.genes = var.genes, do.umap = TRUE)
saveRDS(earlyRetina_all, file="R_objects/earlyRetina_all.rds")

earlyRetina_all = readRDS("R_objects/earlyRetina_all.rds")
earlyRetina_all = set.all.ident(earlyRetina_all, id="age")



earlyRetina_all@data.info$sample = "T"
earlyRetina_all@data.info[grep("Clark", earlyRetina_all@cell.names),"sample"] = "C"
earlyRetina_all@data.info[grep("LoGiudice", earlyRetina_all@cell.names),"sample"] = "L"
earlyRetina_all = set.all.ident(earlyRetina_all, id="sample")

# Rename 
earlyRetina_all@data.info$RGCs_only = "Other_ThisStudy"
earlyRetina_all@data.info[grep("Clark", earlyRetina_all@cell.names),"RGCs_only"] = "Other_Clark"
earlyRetina_all@data.info[which.cells(earlyRetina_all, "ThisStudy_RGC"),"RGCs_only"] = "RGC_ThisStudy"
earlyRetina_all@data.info[which.cells(earlyRetina_all, "Clark_Retinal Ganglion Cells"),"RGCs_only"] = "RGC_Clark"
 DimPlot(earlyRetina_all,group.by = "RGCs_only", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, do.return = TRUE, cols.use = c("gray50","gray70", "blue","cyan"))

p=list()
p[[1]] = DimPlot(earlyRetina_all, reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, cols.use = c("green","blue","gray70"), do.return = TRUE)
p[[2]] = DimPlot(earlyRetina_all, cells.use = which.cells(earlyRetina_all, "C"), reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, cols.use = c("green"), do.return = TRUE)
p[[3]] = DimPlot(earlyRetina_all, cells.use = which.cells(earlyRetina_all, "L"), reduction.use = "UMAP", pt.size = 0.3, cols.use = c("blue"), do.return = TRUE)
p[[4]] = DimPlot(earlyRetina_all, cells.use = which.cells(earlyRetina_all, "T"),  reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, cols.use = c("gray70"), do.return = TRUE)


pdf("Figs/ClarkLG_comparison_UMAP.pdf",w=9,h=7,useDingbats = FALSE)
plot_grid(plotlist = p, nrow=2)
dev.off()

p=list()
p[[1]] = FeaturePlot(earlyRetina_all, "Nefl", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, no.legend = FALSE, do.return = TRUE, cols.use = c("cyan","red"))
p[[1]] = p[[1]][[1]] + theme_classic()
p[[2]] = FeaturePlot(earlyRetina_all, "Fgf15", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, no.legend = FALSE, do.return = TRUE, cols.use = c("cyan","red"))
p[[2]] = p[[2]][[1]] + theme_classic()
p[[3]] = FeaturePlot(earlyRetina_all, "Tfap2b", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, no.legend = FALSE, do.return = TRUE, cols.use = c("cyan","red"))
p[[3]] = p[[3]][[1]] + theme_classic()
p[[4]] = FeaturePlot(earlyRetina_all, "Gngt2", reduction.use = "UMAP", max.cells = 10000, pt.size = 0.3, no.legend = FALSE, do.return = TRUE, cols.use = c("cyan","red"))
p[[4]] = p[[4]][[1]] + theme_classic()

pdf("Figs/ClarkLG_comparison_genes.pdf",w=9,h=7,useDingbats = FALSE)
plot_grid(plotlist = p, nrow=2)
dev.off()

