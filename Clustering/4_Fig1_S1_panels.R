assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# load E13
rgcE13 = readRDS("../Objects/rgcE13_190911_withP56.rds")
# load E14
rgcE14 = readRDS("../Objects/rgcE14_190911_withP56.rds")
# load E16
rgcE16 = readRDS("../Objects/rgcE16_190911_withP56.rds")
# load P0
rgcP0 = readRDS("../Objects/rgcP0_190911_withP56.rds")
# load P5
rgcP5 = readRDS("../Objects/rgcP5_190911_withP56.rds")

# load atlas
rgc_atlas = readRDS("../Objects/rgc_atlas_190911_withP56.rds")
rgc_atlas = set.all.ident(rgc_atlas, id="mnn_final")

objects = c("rgcE13","rgcE14","rgcE16","rgcP0","rgcP5","rgc_atlas")
for (n in c(1:6)){
  eval(parse(text=paste0("object = ", objects[n])))
  fig.name = paste0("Figs/UMAP_", objects[n],".pdf")
  pdf(fig.name, w=6,h=5.5,useDingbats = FALSE)
    DimPlot(object, reduction.use = "UMAP", do.label = TRUE, no.legend = TRUE, pt.size = 0.4)
  dev.off()
}

objects = c("rgcE13","rgcE14","rgcE16","rgcP0","rgcP5","rgc_atlas")
for (n in c(1:6)){
  eval(parse(text=paste0("object = ", objects[n])))
  fig.name = paste0("Figs/UMAP_", objects[n],".pdf")
  pdf(fig.name, w=6,h=5.5,useDingbats = FALSE)
  DimPlot(object, reduction.use = "UMAP", do.label = TRUE, no.legend = TRUE, pt.size = 0.4)
  dev.off()
}

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
# Early Retina
earlyRetina = readRDS("Objects/earlyRetina_E13toP5.rds")
earlyRetina = subsetData(earlyRetina, cells.use = setdiff(earlyRetina@cell.names, which.cells(earlyRetina,"P3")))
genes.use = c("Pou4f1","Rbpms","Syt4","Pou6f2","Slc17a6",
              "Tfap2b","Tfap2a","Onecut2","Ptf1a",  "Otx2","Gngt2","Gnb3","Crx",  "Ctss","C1qa","P2ry12","Fcrls","Ccnd1","Fgf15","Hes5","Hes6",  "Sfrp2", "Mgp","Rgs5","Igfbp7","Col3a1","Bgn")

genes.use = c("Sfrp2","Ccnd2","Nfia","Klf9","Sox9","Nfix","Hes1","Hes5","Id1","Id3","Mki67")

a=DoHeatmap(earlyRetina, use.scaled=TRUE, genes.use=genes.use, group.by = "ident", max.cells.per.ident = 8000,
            ident.order = c("RGC","AC","Cone","Microglia", "P1","P2","P3"), 
            cex.col = 0.01, cex.row=3, col.low = "white", col.mid="white", col.high = "gray40")
pdf("Figs/Heatmap_classes_earlyRetina.pdf", w=10,h=5.5,useDingbats = FALSE)
print(a)
dev.off()

pdf("Figs/UMAP_earlyRetina_classes.pdf", w=7,h=5.5,useDingbats = FALSE)
DimPlot(earlyRetina, reduction.use = "UMAP",  pt.size=0.3)
dev.off()

# Shannon entropy
ShannonH = function(x){y = x[x>0]; y = y/sum(y); return(-sum(y*log2(y)))}
cell_entropy1 = apply(rgc_dev@count.data[,1:50000], 2, function(x) return(ShannonH(x)))
cell_entropy2 = apply(rgc_dev@count.data[,50001:ncol(rgc_dev@count.data)], 2, function(x) return(ShannonH(x)))

rgc_dev@data.info$shannonH = 0
rgc_dev@data.info[names(cell_entropy1),"shannonH"] = cell_entropy1
rgc_dev@data.info[names(cell_entropy2),"shannonH"] = cell_entropy2



pdf("Figs/UMAP_earlyRetina_age.pdf", w=7,h=5.5,useDingbats = FALSE)
DimPlot(earlyRetina, reduction.use = "UMAP",  pt.size=0.3, group.by = "age")
dev.off()

library(Hmisc)
earlyRetina@data.info$broadAnn = drop.levels(earlyRetina@data.info$broadAnn)
pdf("Figs/BroadAnn_dist_age.pdf", w=6,h=5,useDingbats = FALSE)
cluster.dist.by.ident(earlyRetina, ident1.use = "broadAnn", ident2.use = "age", ylab.use = "Fraction", xlab.use = "Age")
dev.off()

pdf("Figs/BroadAnn_dist_reps.pdf", w=8,h=5,useDingbats = FALSE)
cluster.dist.by.ident(earlyRetina, ident1.use = "broadAnn", ident2.use = "orig", ylab.use = "Fraction", xlab.use = "Age")
dev.off()

# All RGCs
rgc_dev = readRDS("../Objects/rgc_dev_190910.rds")
rgc_dev = set.all.ident(rgc_dev, id="age")

genes.use =  c("Rbpms","Pou4f1","Pou4f2", "Slc17a6", "Thy1", "L1cam")
df = as.data.frame(t(as.matrix(rgc_dev@data[genes.use,])))
df$age = rgc_dev@data.info$age

df_mean = df %>% dplyr::group_by(age) %>% dplyr::summarise_all(list( mean))
df_mean = melt(as.data.frame(df_mean))
colnames(df_mean) = c("Age","Gene","Mean")

df_sd = df %>% dplyr::group_by(age) %>% dplyr::summarise_all(list(sd))
df_sd = melt(as.data.frame(df_sd))
colnames(df_sd) = c("Age","Gene","SD")

df_mean$SD = df_sd$SD

pdf("Figs/RGC_markers_with_age.pdf", w=9,h=6, useDingbats = FALSE)
ggplot(df_mean, aes(x=Age, y=Mean)) + geom_bar(stat="identity", fill="darksalmon", color="black") + geom_errorbar(aes(ymin = Mean-SD/2, ymax=Mean+SD/2), width=0.2) + facet_wrap(~ Gene, scales = "free") + theme_classic() + ylab("Expression")
dev.off()

# Diversity
shannonH = rep(0,6); names(shannonH) = levels(rgc_dev@data.info$age)
SimpsonI = shannonH
normShannonH = shannonH
normSimpsonI = SimpsonI
for (n in c(1:6)){
  eval(parse(text=paste0("object = ", objects[n])))
  p = as.numeric(table(object@ident)); p = p/sum(p)
  shannonH[n] = -sum(p*log(p))
  SimpsonI[n] = sum(p^2)
  normShannonH[n] = shannonH[n] / log(length(p))
  normSimpsonI[n] = SimpsonI[n] / (1/length(p))
}

#GiniI[1] = 0.816; GiniI[2] = 0.885
#shannonH[1] = 2.15; shannonH[2] = 2.303

df = data.frame(shannonH = shannonH, SimpsonI = SimpsonI )
df$age = rownames(df)
pdf("Figs/Shannon_Gini.pdf", w=6,h=5, useDingbats = FALSE)
ggplot(df, aes(x=age)) + geom_point(aes(y=shannonH, color="Shannon"), size=2) + geom_line(aes(y=shannonH, group=1, color="Shannon"), size = 1.5) + geom_point(aes(y=40*SimpsonI,  color="Simpson"), size=2) + geom_line(aes(y=SimpsonI*40, group=1, color="Simpson"), size=1.5) +  scale_y_continuous(sec.axis = sec_axis(~./40, name = "Simpson Index"), limits = c(0,5)) + theme_classic() + scale_color_manual(values = c("blue", "red")) + theme(legend.position = "top") + labs(x = "Age", y = "Shannon Entropy")
dev.off()

rgc_atlas@data.info$all = 1
ids = c(rep("m_merge",4), rep("mnn_final",2))
# Correlations 
for (n in c(1:6)){
  eval(parse(text=paste0("object = ", objects[n])))
  
  object = set.all.ident(object, id="all")
  var.genes500 = NB.var.genes(object, x.low.cutoff = 0.01, x.high.cutoff = 30, return.top.genes = TRUE, topn = 500 )
 object = set.all.ident(object, id=ids[n]) 
  A = Avg.by.ident(object,features.use = names(var.genes500), return.scaled = FALSE )
  C = cor(A)
  vals = C[upper.tri(C)]
  
  if (n == 1){
    df = data.frame(cor_vals = vals, age = levels(rgc_dev@ident)[n])
  } else {
    df = rbind(df, data.frame(cor_vals = vals, age = levels(rgc_dev@ident)[n]))
  }
}

df_mean = as.data.frame(df %>% dplyr::group_by(age) %>% summarize_all(median))


# Clusters vs. cells 
clusters = c(10,16, 19, 27, 38, 45)
cells = c(5973, 17100, 13046, 18148, 17386, 35699 )

