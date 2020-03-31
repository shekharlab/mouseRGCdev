assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Read Count matrices
Count_mat_E13toP0 = readRDS("../ConsolidatedCounts/Counts_E13toP0.rds")
Count_mat_p5 = readRDS("../../ConsolidatedCounts/Counts_P5.rds")
Count_mat_p5 = Count_mat_p5$All

# Filter Count matrix
genes.common = intersect(rownames(Count_mat_E13toP0), rownames(Count_mat_p5))
Count.mat = cbind(Count_mat_E13toP0[genes.common,], Count_mat_p5[genes.common,])
ngenes = colSums(Count.mat > 0)
Count.mat = Count.mat[,ngenes > 700]

# Create S4 object
earlyRetina=scR(count.data=Count.mat,ident.fxn=getStat1)
earlyRetina=setup(earlyRetina,project="earlyRetina",min.cells = 10,min.genes = 700,is.expr=0, threshold.quantile = 0.999)

# Total 98,452 cells
earlyRetina@data = Matrix(earlyRetina@data, sparse=TRUE)

earlyRetina@data.info$age = "P5"
earlyRetina@data.info[grep("^E13", earlyRetina@cell.names),"age"] = "E13"
earlyRetina@data.info[grep("^E14", earlyRetina@cell.names),"age"] = "E14"
earlyRetina@data.info[grep("^E16", earlyRetina@cell.names),"age"] = "E16"
earlyRetina@data.info[grep("^P0", earlyRetina@cell.names),"age"] = "P0"

earlyRetina@data.info$age = factor(earlyRetina@data.info$age, levels = c("E13","E14", "E16", "P0", "P5"))
earlyRetina = set.all.ident(earlyRetina, id="age")
VlnPlot(earlyRetina, c("nGene","Rbpms"),nCol=1)
earlyRetina = set.all.ident(earlyRetina, id="orig")
VlnPlot(earlyRetina, c("nGene","Rbpms"), nCol=1, x.lab.rot=TRUE)

# Ligerize
earlyRetina@data.info$all = 1
earlyRetina = set.all.ident(earlyRetina, id="all")
var.genes = NB.var.genes(earlyRetina,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0001, rem.mt.rp = FALSE, do.text = FALSE)
earlyRetina@var.genes = var.genes

earlyRetina = Ligerize(earlyRetina, batch_id = "age", var.genes = var.genes, do.umap = TRUE)

# Define cell classes
earlyRetina = buildClusterTree(earlyRetina, genes.use = earlyRetina@var.genes, linkage.method = "complete", dist.fun = "euclidean", do.scale=FALSE, regress.use = FALSE)

pdf("Figs/earlyRetina_all_dendro_dotplot.pdf",w=13,h=7, useDingbats=FALSE)
plot_dendro_withdotplot(earlyRetina,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25"), 
                        max.val.exp = 5, alpha.use=0.7)
dev.off()

pdf("Figs/earlyRetina_dendro.pdf",w=13,h=7, useDingbats=FALSE)
plotClusterTree(earlyRetina)
dev.off()

plot_dendro_with_sample_dist(earlyRetina, id.use = "age")

# Define classes
precursor1= setdiff(getLeftDecendants(earlyRetina@cluster.tree[[1]],49),c(22,30))
precursor2 = c(22,30)
rgc.clusters = c(getLeftDecendants(earlyRetina@cluster.tree[[1]],51), getLeftDecendants(earlyRetina@cluster.tree[[1]],58))
mic.cluster = 24
cone.prec = 28
ac.clusters = c(9,37)
precursor3 = getRightDecendants(earlyRetina@cluster.tree[[1]],58)

all_clusters = c(precursor1, precursor2, precursor3,
                 rgc.clusters, mic.cluster, cone.prec, ac.clusters)

# P1 = RPCs, P2 = Anterior Segment Cells, P3 = Doublets
earlyRetina@data.info$broadAnn = "RGC"
earlyRetina@data.info[which.cells(earlyRetina, precursor1),"broadAnn"] = "P1"
earlyRetina@data.info[which.cells(earlyRetina, precursor2),"broadAnn"] = "P2"
earlyRetina@data.info[which.cells(earlyRetina, precursor3),"broadAnn"] = "P3"
earlyRetina@data.info[which.cells(earlyRetina, ac.clusters),"broadAnn"] = "AC"
earlyRetina@data.info[which.cells(earlyRetina, mic.cluster),"broadAnn"] = "Microglia"
earlyRetina@data.info[which.cells(earlyRetina, cone.prec),"broadAnn"] = "Cone"

earlyRetina@data.info$broadAnn = factor(earlyRetina@data.info$broadAnn,
                                       levels= c("Microglia","Cone","AC","P1","P2","P3","RGC"))


# Extract RGCs
earlyRetina = set.all.ident(earlyRetina, id="broadAnn")
earlyRGC = subsetData(earlyRetina, cells.use = which.cells(earlyRetina,"RGC"))
earlyRGC = set.all.ident(earlyRGC, id="age")
rgcE13 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "E13"))
rgcE14 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "E14"))
rgcE16 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "E16"))
rgcP0 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "P0"))
rgcP5 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "P5"))

# Save Objects
saveRDS(earlyRetina, file="../R_objects/earlyRetina_E13toP5.rds")
saveRDS(earlyRGC, file="../R_objects/earlyRGC_E13toP5.rds")
saveRDS(rgcE13, file="../R_objects/rgcE13_190729.rds")
saveRDS(rgcE14, file="../R_objects/rgcE14_190729.rds")
saveRDS(rgcE16, file="../R_objects/rgcE16_190729.rds")
saveRDS(rgcP0, file="../R_objects/rgcP0_190729.rds")
saveRDS(rgcP5, file="../R_objects/rgcP5_190729.rds")


#########################################################################
###### Fig 1,2 and S1 panels
#######################################################################
library(Hmisc)

### Figure 1, S1
earlyRetina = set.all.ident(earlyRetina, id = "age")
dir.create("Figs")

# Fig 1C, D
pdf("Figs/UMAP_earlyRGCs_age.pdf",w=8,h=6, useDingbats=FALSE)
DimPlot(earlyRetina, cells.use = sample(earlyRetina@cell.names, 20000), reduction.use = "UMAP", pt.size = 0.3)
dev.off()
earlyRetina = set.all.ident(earlyRetina, id="broadAnn")
pdf("Figs/UMAP_earlyRGCs_broadAnn.pdf",w=8,h=6, useDingbats=FALSE)
DimPlot(earlyRetina, cells.use = sample(earlyRetina@cell.names, 20000), reduction.use = "UMAP", pt.size = 0.3, do.label=TRUE)
dev.off()


# Identify cell classes and remove contaminants
earlyRetina = find.all.markers(earlyRetina,test.use="MAST", max.cells.per.ident = 2000)
earlyRetina = buildClusterTree(earlyRetina, genes.use=earlyRetina@var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "correlation")


p = Perc.pos.by.ident(earlyRetina, c("Ctss","C1qa","P2ry12","Gngt2","Otx2","Gnb3",
                                     "Tfap2b","Tfap2a","Onecut2","Ccnd1","Fgf15","Hes5","Sfrp2","Mgp","Rgs5","Igfbp7","Col3a1","Pou4f1","Rbpms","Sncg", "Slc17a6","Atoh7"),
                      ident.use = c("Microglia", "Cone", "AC", "P1","P2", "RGC"), x.lab.rot=TRUE, norm.exp = 1,return.plot=TRUE, alpha.use = 0.7)

# Figure X (Not used)
pdf("Figs/CellClasses_earlyRetina.pdf",w=4,h=7,useDingbats = FALSE)
p
dev.off()

# Markers between progenitors
a=find.markers(earlyRetina,"P1","P2", test.use="MAST", max.cells.per.ident = 1000)
p = Perc.pos.by.ident(earlyRetina, c(rownames(head(subset(a,diff>0),10)), rownames(tail(subset(a,diff<0),10))),ident.use = c("P1","P2"), x.lab.rot=TRUE, norm.exp = 1,return.plot=TRUE, alpha.use = 0.7)

# Not used
pdf("Figs/P1P2_earlyRetina.pdf",w=3,h=8,useDingbats = FALSE)
p
dev.off()


# Figure S1B
earlyRetina@data.info$age_sort = as.character(earlyRetina@data.info$orig)
for (i in c(1:16)){
  earlyRetina@data.info$age_sort = gsub(paste0("RGCS",i,"$"),"",earlyRetina@data.info$age_sort)
  
}
pdf("Figs/earlyRetina_broadAnn_dist_groups.pdf",w=6,h=5, useDingbats = FALSE)
cluster.dist.by.ident(earlyRetina, ident1.use = "broadAnn", ident2.use = "age_sort")
dev.off()


# Heatmap (Fig 1E)
genes.use = c("Snap25","Nefl", "Pou4f1","Rbpms","Syt4","Pou6f2","Slc17a6",
              "Tfap2b","Tfap2a","Onecut2","Ptf1a",  "Otx2","Gngt2","Gnb3","Crx",  "Ctss","C1qa","P2ry12","Fcrls","Ccnd1","Fgf15","Hes5","Hes6",  "Sfrp2", "Mgp","Rgs5","Igfbp7","Col3a1","Bgn","Aldh1a1","Celf4")
a=DoHeatmap(earlyRetina, use.scaled=TRUE, genes.use=genes.use, group.by = "ident", max.cells.per.ident = 2000,
            ident.order = c("RGC","AC","Cone","Microglia", "P1","P2","P3"), 
            cex.col = 0.01, cex.row=3, col.low = "white", col.mid="white", col.high = "red")
pdf("Figs/Heatmap_classes_earlyRetina.pdf", w=10,h=5.5,useDingbats = FALSE)
print(a)
dev.off()

# Figure S1D
objects2 = c("rgcE14","rgcE16","rgcP0")
p=list()
for (n in c(1:3)){
  eval(parse(text=paste0("object = ", objects2[n])))
  object@data.info$age_sort = as.character(object@data.info$orig)
  for (i in c(1:16)){
    object@data.info$age_sort = gsub(paste0("RGCS",i,"$"),"",object@data.info$age_sort)
  }
  
  object = set.all.ident(object, id="age_sort")
  p[[n]] = DimPlot(object, reduction.use = "UMAP", do.label = FALSE, pt.size = 0.4, do.return = TRUE, cols.use = c("gray20","gray80"))
  
}

pdf("Figs/UMAP_L1cam_Cd90.pdf", w=6,h=12,useDingbats = FALSE)
plot_grid(plotlist=p, nrow=3)
dev.off()

