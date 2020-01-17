assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Read Count matrices
Count_mat_E13toP0 = readRDS("../../RawData/Counts_E13toP0.rds")
Count_mat_p5 = readRDS("../../RawData/Counts_P5.rds")
Count_mat_p5 = Count_mat_p5$All

# Filter Count matrix
genes.common = intersect(rownames(Count_mat_E13toP0), rownames(Count_mat_p5))
Count.mat = cbind(Count_mat_E13toP0[genes.common,], Count_mat_p5[genes.common,])
ngenes = colSums(Count.mat > 0)
Count.mat = Count.mat[,ngenes > 700]

# Create object
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

dir.create("../Objects")
saveRDS(earlyRetina, file="../Objects/earlyRetina_E13toP5.rds")
earlyRetina = readRDS("../Objects/earlyRetina_E13toP5.rds")

# Ligerize
earlyRetina@data.info$all = 1
earlyRetina = set.all.ident(earlyRetina, id="all")
var.genes = NB.var.genes(earlyRetina,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0001, rem.mt.rp = FALSE, do.text = FALSE)
earlyRetina@var.genes = var.genes

earlyRetina = Ligerize(earlyRetina, batch_id = "age", var.genes = var.genes, do.umap = TRUE)
saveRDS(earlyRetina, file="../Objects/earlyRetina_E13toP5.rds")

dir.create("Figs")
pdf("Figs/UMAP_earlyRGCs_age.pdf",w=10,h=8, useDingbats=FALSE)
DimPlot(earlyRetina, reduction.use = "UMAP", pt.size = 0.3)
dev.off()
earlyRetina = set.all.ident(earlyRetina, id="m")
pdf("Figs/UMAP_earlyRGCs_cluster.pdf",w=11,h=8, useDingbats=FALSE)
DimPlot(earlyRetina, reduction.use = "UMAP", pt.size = 0.3, do.label = TRUE)
dev.off()

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

earlyRetina@data.info$broadAnn = "RGC"
earlyRetina@data.info[which.cells(earlyRetina, precursor1),"broadAnn"] = "P1"
earlyRetina@data.info[which.cells(earlyRetina, precursor2),"broadAnn"] = "P2"
earlyRetina@data.info[which.cells(earlyRetina, precursor3),"broadAnn"] = "P3"
earlyRetina@data.info[which.cells(earlyRetina, ac.clusters),"broadAnn"] = "AC"
earlyRetina@data.info[which.cells(earlyRetina, mic.cluster),"broadAnn"] = "Microglia"
earlyRetina@data.info[which.cells(earlyRetina, cone.prec),"broadAnn"] = "Cone"

earlyRetina@data.info$broadAnn = factor(earlyRetina@data.info$broadAnn,
                                       levels= c("Microglia","Cone","AC","P1","P2","P3","RGC"))



pdf("Figs/earlyRetina_broadAnn_dist_age.pdf",w=6,h=4, useDingbats = FALSE)
cluster.dist.by.ident(earlyRetina, ident2.use = "age", ident1.use = "broadAnn")
dev.off()

earlyRetina = set.all.ident(earlyRetina, id="broadAnn")

# Monocle
gene_annotation = data.frame(gene_short_name = rownames(earlyRetina@count.data));
rownames(gene_annotation) = gene_annotation$gene_short_name
cds <- new_cell_data_set(earlyRetina@count.data,
                         cell_metadata = earlyRetina@data.info,
                         gene_metadata = gene_annotation)
reducedDims(cds)$UMAP = earlyRetina@dr$UMAP@cell.embeddings
marker_test_res = top_markers(cds, group_cells_by="broadAnn", reference_cells=1000, cores=8)
top_specific_markers = marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(4, pseudo_R2)

top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="broadAnn",
                    ordering_type="maximal_on_diag",
                    max.size=3)

earlyRetina = find.all.markers(earlyRetina,test.use="MAST", max.cells.per.ident = 2000)

earlyRetina = buildClusterTree(earlyRetina, genes.use=earlyRetina@var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "correlation")

p = Perc.pos.by.ident(earlyRetina, c("Ctss","C1qa","P2ry12","Gngt2","Otx2","Gnb3",
                                     "Tfap2b","Tfap2a","Onecut2","Ccnd1","Fgf15","Hes5","Sfrp2","Mgp","Rgs5","Igfbp7","Col3a1","Pou4f1","Rbpms","Sncg", "Slc17a6","Atoh7"),
                      ident.use = c("Microglia", "Cone", "AC", "P1","P2", "RGC"), x.lab.rot=TRUE, norm.exp = 1,return.plot=TRUE, alpha.use = 0.7)

pdf("Figs/CellClasses_earlyRetina.pdf",w=4,h=7,useDingbats = FALSE)
p
dev.off()

# Markers between progenitors
a=find.markers(earlyRetina,"P1","P2", test.use="MAST", max.cells.per.ident = 1000)
p = Perc.pos.by.ident(earlyRetina, c(rownames(head(subset(a,diff>0),10)), rownames(tail(subset(a,diff<0),10))),ident.use = c("P1","P2"), x.lab.rot=TRUE, norm.exp = 1,return.plot=TRUE, alpha.use = 0.7)

pdf("Figs/P1P2_earlyRetina.pdf",w=3,h=8,useDingbats = FALSE)
p
dev.off()

VlnPlot(earlyRetina,"nGene", x.lab.rot=TRUE)
saveRDS(earlyRetina, file="../Objects/earlyRetina_E13toP5.rds")
earlyRetina = readRDS("../Objects/earlyRetina_E13toP5.rds")

earlyRetina = subsetData(earlyRetina, cells.use = setdiff(earlyRetina@cell.names,
                                                          which.cells(earlyRetina,"P3")))

pdf("Figs/UMAP_earlyRGCs_broadAnn.pdf",w=10,h=8, useDingbats=FALSE)
DimPlot(earlyRetina, reduction.use = "UMAP", pt.size = 0.3, do.label=TRUE)
dev.off()

earlyRetina@data.info$broadAnn = drop.levels(earlyRetina@data.info$broadAnn)
library(Hmisc)
pdf("Figs/earlyRetina_broadAnn_dist.pdf",w=12,h=6, useDingbats = FALSE)
cluster.dist.by.ident(earlyRetina, ident1.use = "broadAnn")
dev.off()

# Look at P2
P2 = subsetData(earlyRetina, cells.use = which.cells(earlyRetina,"P2"))
P2 = set.all.ident(P2, id="age")
a = find.markers(P2, "E13","E14", test.use="MAST", max.cells.per.ident = 1000)

# Look at P1
P1 = subsetData(earlyRetina, cells.use = which.cells(earlyRetina,"P1"))
P1 = set.all.ident(P1, id="age")
a = find.markers(P1, "E14","P0", test.use="MAST", max.cells.per.ident = 1000)

# Extract RGCs
earlyRetina = set.all.ident(earlyRetina, id="broadAnn")
earlyRGC = subsetData(earlyRetina, cells.use = which.cells(earlyRetina,"RGC"))
saveRDS(earlyRGC, file="../Objects/earlyRGC_E13toP5.rds")
earlyRGC = readRDS("../Objects/earlyRGC_E13toP5.rds")
earlyRGC = set.all.ident(earlyRGC, id="age")
rgcE13 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "E13"))
rgcE14 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "E14"))
rgcE16 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "E16"))
rgcP0 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "P0"))
rgcP5 = subsetData(earlyRGC, cells.use = which.cells(earlyRGC, "P5"))

saveRDS(rgcE13, file="../Objects/rgcE13_190729.rds")
saveRDS(rgcE14, file="../Objects/rgcE14_190729.rds")
saveRDS(rgcE16, file="../Objects/rgcE16_190729.rds")
saveRDS(rgcP0, file="../Objects/rgcP0_190729.rds")
saveRDS(rgcP5, file="../Objects/rgcP5_190729.rds")
