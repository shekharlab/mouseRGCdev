assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Load full RGC object
rgc_dev = readRDS("../Developing_RGCs/Objects/rgc_dev_190910.rds")

# Load individual objects
# load E13
rgcE13 = readRDS("../Developing_RGCs/Objects/rgcE13_190911_withP56.rds")

# load E14
rgcE14 = readRDS("../Developing_RGCs/Objects/rgcE14_190911_withP56.rds")

# load E16
rgcE16 = readRDS("../Developing_RGCs/Objects/rgcE16_190911_withP56.rds")

# load P0
rgcP0 = readRDS("../Developing_RGCs/Objects/rgcP0_190911_withP56.rds")

# load P5
rgcP5 = readRDS("../Developing_RGCs/Objects/rgcP5_190911_withP56.rds")

# Load atlas
rgc_atlas=readRDS("../Developing_RGCs/Objects/rgc_atlas_190911_withP56.rds")




objects = c("rgcE13","rgcE14","rgcE16","rgcP0","rgcP5","rgc_atlas")
traj_files = c("e13_traj", "e14_traj","e16_traj","p0_traj","p5_traj","p56_traj")
ages = c("E13", "E14","E16","P0","P5","P56")

rgc_atlas = set.all.ident(rgc_atlas, id="new_types")
# Compute fate clusters at each time point
for (k in c(1:5)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  
  # Use trajectory values
  data.use = object@data.info[,paste0(levels(rgc_atlas@ident),"_fate")]
  Adj = getAdjMatrix(data.use,nn=40,edge.weights=FALSE,do.jaccard=TRUE,do.sparse=TRUE,full.eval=FALSE)
  g <- graph.adjacency(Adj, weighted=TRUE, mode="undirected")
  graph.out = cluster_louvain(g)
  object@data.info$fate_clusters = factor(graph.out$membership)
  eval(parse(text = paste0(objects[k], "= object")))
}


saveRDS(rgcE13, file="../Developing_RGCs/Objects/rgcE13_190911_withP56.rds")
saveRDS(rgcE14, file="../Developing_RGCs/Objects/rgcE14_190911_withP56.rds")
saveRDS(rgcE16, file="../Developing_RGCs/Objects/rgcE16_190911_withP56.rds")
saveRDS(rgcP0, file="../Developing_RGCs/Objects/rgcP0_190911_withP56.rds")
saveRDS(rgcP5, file="../Developing_RGCs/Objects/rgcP5_190911_withP56.rds")

# Write Cell Sets for OT
# Write cell sets
objects2 = c("rgc_atlas","rgcP5", "rgcP0", "rgcE16", "rgcE14","rgcE13")
ages2 = c("P56","P5","P0","E16","E14","E13")
sink("OT_uncorr_vargenes_eps005/cell_sets_fate_clusters.gmt")

#library(xgboost)
for (n in c(1:6)){
  
  eval(parse(text=paste0("object = ", objects2[n])))
  if (n > 1) object = set.all.ident(object, id="fate_clusters")
  for (i in levels(object@ident)){
    cells.use = which.cells(object, i)
    if (n != 1){
      k = paste0("C",i)
    } else {
      k = i
    }
    set_name = paste0(ages2[n], "_", k)
    cells_set = paste0(cells.use, collapse = "\t")
    
    str_to_write = paste0(set_name,"\t NA \t", cells_set)
    cat(str_to_write)
    cat("\n")
  }
}
sink()


# Plot Confusion Matrix between fate clusters and transcriptional clusters
plist = list()
id.list = c("m_merge","m_merge","m_merge","m_merge","mnn_final")
library(infotheo)
for (k in c(1:5)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  A = table(object@data.info$fate_clusters, object@data.info[,id.list[k]])
  test <- heatmap.2(A, hclustfun = function(x) hclust(x,method="single"))
  B=A[test$rowInd, test$colInd]
  row.max = apply(B,1,which.max)
  row.ord = names(sort(row.max))
  NCE = condentropy(object@data.info[,id.list[k]],object@data.info$fate_clusters) / entropy(object@data.info[,id.list[k]])
  
  if (k > 2) NCE = NCE*0.75
  print(NCE)
plist[[k]]=plotConfusionMatrix(B[row.ord,], col.low = "white", col.high = "darkblue", ylab.use = paste0(ages[k], " fate clusters"), xlab.use = paste0(ages[k]," transcriptional clusters"), plot.return = TRUE) + ggtitle(paste0("NCE = ", round(NCE,2)))
}

pdf("Figs_3_4/ConfusionMatrices_FateVsTranscriptionalClusters.pdf",w=35,h=8, useDingbats = FALSE)
plot_grid(plotlist = plist[1:5], nrow=1, rel_widths = c(1,1.2,1.4,1.6,2))
dev.off()

pdf("Figs_3_4/ConfusionMatrices_FateVsTranscriptionalClusters_col.pdf",w=8,h=32, useDingbats = FALSE)
plot_grid(plotlist = plist[1:5], ncol=1)
dev.off()


# UMAP
plist = list()
id.list = c("m_merge","m_merge","m_merge","m_merge","mnn_final")
library(infotheo)
for (k in c(1:5)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  plist[[k]] = DimPlot(object, reduction.use="UMAP", pt.size=0.3, do.return = TRUE, group.by = "fate_clusters")
}

pdf("Figs_3_4/UMAP_fate_clusters_col.pdf",w=7,h=32, useDingbats = FALSE)
plot_grid(plotlist = plist[1:5], ncol=1)
dev.off()

# Fate bias for each cluster
plist = list()
for (k in c(1:5)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  data.use = object@data.info[,c(paste0(levels(rgc_atlas@ident),"_fate"), "fate_clusters")] 
  data.avg = data.use %>% dplyr::group_by(fate_clusters) %>% dplyr::summarise_all(funs(mean))
  data.avg = as.data.frame(data.avg)
  data.avg = data.avg[,-1]
  data.avg = scale(data.avg, center=FALSE, scale = colSums(data.avg))
  rownames(data.avg) = levels(object@data.info$fate_clusters);
  colnames(data.avg) = c(1:45)
  
  test <- heatmap.2(data.avg, hclustfun = function(x) hclust(x,method="single"))
  B=data.avg[test$rowInd, test$colInd]
  row.max = apply(B,1,which.max)
  row.ord = names(sort(row.max))
  plist[[k]] = DataHeatmap2(t(B[row.ord,]), vals.range = c(1e-5,0.5), 
               col.low = "white", col.mid = "white", col.high = "darkgreen", return.plot = TRUE)
  
}

pdf("Figs_3_4/Type_bias_fate_clusters.pdf",w=35,h=7, useDingbats = FALSE)
plot_grid(plotlist = plist[1:5], nrow=1)
dev.off()


# Fate bias for each transcriptional cluster
plist = list()
for (k in c(1:5)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  data.use = object@data.info[,c(paste0(levels(rgc_atlas@ident),"_fate"), id.list[k])] 
  colnames(data.use)[ncol(data.use)] = "clusters"
  data.avg = data.use %>% dplyr::group_by(clusters) %>% dplyr::summarise_all(funs(mean))
  data.avg = as.data.frame(data.avg)
  data.avg = data.avg[,-1]
  data.avg = scale(data.avg, center=FALSE, scale = colSums(data.avg))
  rownames(data.avg) = levels(object@data.info[,id.list[k]]);
  colnames(data.avg) = c(1:45)
  
  test <- heatmap.2(data.avg, hclustfun = function(x) hclust(x,method="single"))
  B=data.avg[test$rowInd, test$colInd]
  row.max = apply(B,1,which.max)
  row.ord = names(sort(row.max))
  plist[[k]] = DataHeatmap2(t(B[row.ord,]), vals.range = c(1e-5,0.5), 
                            col.low = "white", col.mid = "white", col.high = "darkgreen", return.plot = TRUE)
  
}

pdf("Figs_3_4/Type_bias_trans_clusters.pdf",w=35,h=7, useDingbats = FALSE)
plot_grid(plotlist = plist[1:5], nrow=1)
dev.off()

library(igraph)
library(qgraph)
# Node colors
cols.use = rep("white",45)
cols.use[c(22,33,40,31,43)] = "blue" # ip RGCs
cols.use[c(42,43,41,45)] = "purple" # Alpha RGCs
cols.use[c(1,2,23,3,4,6,21,30)] = "chartreuse" # T5-RGCs
cols.use[c(38,3,4,28,32)] = "yellow2" # Foxp2
cols.use[c(17,5,9,21,42)] = "orange" # Tbr1
cols.use[c(14,16,36,12)] = "salmon" # Cartpt
cols.use[c(26,19,20,29,35,12,25,39)] = "red" # neurod2


# Correlation between type trajectories
clust.order = c(15,17,5,9,21,13,11,38,28,4,3,6,1,2,32,26,19,20,29,35,12,25,39,44,33,40,22,31,7,8,
                43,42,45,41,27,18,37,36,14,16,24,10,34,30,23)
for (k in c(1:5)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  C = cor(object@data.info[,paste0(levels(rgc_atlas@ident),"_fate")])
  diag(C) = 0
  colnames(C) = c(1:45)
  rownames(C) = c(1:45)
  g <- graph.adjacency(C, weighted=TRUE, mode="undirected")
  minC <- rep(-Inf, vcount(g))
  maxC <- rep(Inf, vcount(g))
 
  g <- graph.adjacency(C, weighted=TRUE, mode="undirected")
  
  C2=C
  C2[C2 < 0.2] = 0
  #C2[C2 < 0.2 & C2 > 0] = 0
  #C2[C2 > -0.2 & C2 < 0] = 0
  g2 <- graph.adjacency(C2, weighted=TRUE, mode="undirected")
  
  if (k == 1){
    set.seed(1)
    #co <- layout_with_fr(g2, maxx = 10000, minx = -10000, maxy = 10000, miny = -10000)
    e <- get.edgelist(g2)
    e2 = matrix(0, nrow=nrow(e), ncol=2)
    e2[,1] = as.numeric(e[,1]); e2[,2] = as.numeric(e[,2])
    set.seed(10)
    co <- qgraph.layout.fruchtermanreingold(e2,vcount=vcount(g2), weights = E(g2)$weight)
  }
  
  
  
 
  C2=C
  C2[C2 < 0.1] = 0
  g3 <- graph.adjacency(C2, weighted=TRUE, mode="undirected")
  widths = 5*E(g3)$weight
  widths[widths > 1] = 1
  
  # Edge colors (color transcriptionally adjacent types in the P56 dendrogram)
  e <- get.edgelist(g3)
  e2 = matrix(0, nrow=nrow(e), ncol=2)
  e2[,1] = as.numeric(e[,1]); e2[,2] = as.numeric(e[,2])
  edge.colors = rep("gray",nrow(e))
  for (ddord in c(2:45)){
    n1 = min(clust.order[ddord-1],clust.order[ddord]); n2 = max(clust.order[ddord-1],clust.order[ddord])
    ind = which(e2[,1] == n1 & e2[,2] == n2)
    if (length(ind) != 0) edge.colors[ind] = "blue"
  }
  
  fig.name = paste0("Figs_3_4/FateCorrelationTypes_age",ages[k],".pdf")
  pdf(fig.name, w=6,h=6, useDingbats = FALSE)
  plot(g3, layout=co, vertex.size=400,vertex.label.cex=0.8,
       vertex.label=c(1:45), rescale=FALSE, edge.width = widths,
       xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=2,
       vertex.label.color="black", asp=1, edge.curved=0.5, vertex.color = cols.use, 
       edge.color=edge.colors)
  dev.off()
  
}

edges3 <- data.frame(Parent = "MP", 
                     Identity = levels(rgc_atlas@ident))

# Muller plot clonal abundance
library(ggmuller)
generation = c(0,1,2,3,4,5)
for (k in c(1:6)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  F = object@data.info[,paste0(levels(rgc_atlas@ident),"_fate")]
  nums = apply(F,2,function(x) sum(x > 0.5))
  nums = c(nrow(F) - sum(nums),nums)
  names(nums) = c("MP", levels(rgc_atlas@ident))
  if (k==1){
    df = data.frame(Generation=generation[k],Identity = names(nums), Population = nums)
  } else {
    df = rbind(df,data.frame(Generation=generation[k],Identity = names(nums), Population = nums))
  }
}

df$Identity = factor(df$Identity, levels = c(levels(rgc_atlas@ident),"MP"))
Muller_df <- get_Muller_df(edges3, df)
Muller_df$Identity = factor(Muller_df$Identity, levels = c(levels(rgc_atlas@ident),"MP"))

pdf("Figs_3_4/MullerPlot.pdf",w=12,h=7,useDingbats = FALSE)
Muller_plot(Muller_df, add_legend = TRUE, xlab = "Time", ylab = "Proportion") + theme_classic()
dev.off()

# Time of emergence
gen_time = c()
for (type in levels(rgc_atlas@ident)){
  df = subset(Muller_df, Identity == type)
  max_freq = 2*max(df$Frequency)
  df_sub = subset(df, Frequency > 0.1*max_freq)
  temp = min(intersect(df_sub$Generation, generation))
  gen_time = c(gen_time, temp)
}

freq_rgc = as.numeric(table(rgc_atlas@ident)/sum(table(rgc_atlas@ident)))

df_plot = data.frame(freq_rgc = freq_rgc, type = levels(rgc_atlas@ident), age = ages[gen_time+1], age=gen_time, subclass = "Other", stringsAsFactors = FALSE)
df_plot[c(22,33,40,31,43),"subclass"] = "ipRGCs"
df_plot[c(42,43,41,45),"subclass"] = "AlphaRGCs"
df_plot[c(1,2,23,3,4,6,21,30),"subclass"] = "T5_RGCs"
df_plot[c(38,3,4,28,32),"subclass"] = "F_RGCs"
df_plot[c(17,5,9,21,42),"subclass"] = "Tbr1_RGCs"
df_plot[c(14,16,36,12),"subclass"] = "Cartpt_RGCs"
df_plot[c(26,19,20,29,35,12,25,39),"subclass"] = "N_RGCs"
df_plot$subclass = factor(df_plot$subclass, levels = c("ipRGCs","AlphaRGCs","T5_RGCs","F_RGCs",
                                                       "Cartpt_RGCs","Tbr1_RGCs","N_RGCs","Other"))



df_plot$fill = cols.use

pdf("Figs_3_4/TypeEmergence_dotplot.pdf",w=7,h=5, useDingbats = FALSE)
ggplot(df_plot, aes(x = factor(age), y = freq_rgc)) +
  geom_dotplot(binaxis = "y", stackdir = "center", stackgroups=TRUE, binpositions="all", aes(fill = factor(subclass)))  + scale_fill_manual(values = c("blue","purple","chartreuse", "yellow2","orange","salmon","red","white")) + scale_y_log10() + xlab("Age") + ylab("Adult Freq") + theme_bw()
dev.off()

