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


# Compute differentially expressed genes across fate clusters

objects = c("rgcE13","rgcE14","rgcE16","rgcP0","rgcP5","rgc_atlas")
ages = c("E13", "E14","E16","P0","P5","P56")
id.list = c("m_merge","m_merge","m_merge","m_merge","mnn_final")

# Compute DE for each fate cluster
for (k in c(1:5)){
  print(k)
  eval(parse(text = paste0("object = ", objects[k])))
  object = set.all.ident(object, id="fate_clusters")
  object = find.all.markers(object, test.use = "MAST", max.cells.per.ident = 1000,add.to.slot = TRUE  )
  names(object@de.list)[length(object@de.list)] = "FateClusterDE_MAST"
  eval(parse(text = paste0(objects[k], "= object")))
}

rgc_atlas = set.all.ident(rgc_atlas, id="mnn_final")
rgc_atlas = find.all.markers(rgc_atlas, test.use = "MAST", max.cells.per.ident = 1000,add.to.slot = TRUE)
names(rgc_atlas@de.list)[length(rgc_atlas@de.list)] = "FateClusterDE_MAST"

# Save objects
saveRDS(rgcE13, file="../Developing_RGCs/Objects/rgcE13_190911_withP56.rds")
saveRDS(rgcE14, file="../Developing_RGCs/Objects/rgcE14_190911_withP56.rds")
saveRDS(rgcE16, file="../Developing_RGCs/Objects/rgcE16_190911_withP56.rds")
saveRDS(rgcP0, file="../Developing_RGCs/Objects/rgcP0_190911_withP56.rds")
saveRDS(rgcP5, file="../Developing_RGCs/Objects/rgcP5_190911_withP56.rds")
saveRDS(rgc_atlas, file="../Developing_RGCs/Objects/rgc_atlas_190911_withP56.rds")


N_vec = c(2,2,2,2,2)
plist1 = list()
plist2 = list()
fate_markers = list()
for (k in c(1:5)){
  print(k)
  N = N_vec[k]
  eval(parse(text = paste0("object = ", objects[k])))
  object = set.all.ident(object, id="fate_clusters")
  df = object@de.list$FateClusterDE_MAST
  df$DEratio = df[,3] * df[,7] / ((df[,4]+0.15) * (df[,8]+0.2))
  ExpMat0 = Perc.pos.by.ident(object, features.use = df$gene, do.plot = FALSE)
  E0 = ExpMat0$ExpMat * ExpMat0$PercMat
  
  markers = c()
  for (clust in 1:length(levels(object@ident))){
    df_sub = subset(df, cluster == clust)
    df_sub = df_sub[order(-df_sub$DEratio),]
    df_sub = subset(df_sub, !(gene %in% markers))
    if (nrow(df_sub) > 0){
      markers = unique(c(markers,df_sub$gene[1:min(N, nrow(df_sub))]))
      # Find the most specific markers ranked by percentage
    }
    if (length(setdiff(df_sub$gene, markers)) > 1){
      ratio_enrich = apply(E0[setdiff(df_sub$gene, markers), ],1,function(x) min(x[clust] / (x[-clust] + 0.05) ))
      ratio_enrich = sort(ratio_enrich, decreasing=TRUE)
      markers = unique(c(markers,names(ratio_enrich)[1:min(N, length(ratio_enrich))]))
    }
    
   
  }
  markers = markers[!is.na(markers)]
  fate_markers[[ages[k]]] = markers
  # Heatmap for fate clusters
  object = set.all.ident(object, id="fate_clusters")
    # Compute Expression matrix
    ExpMat1 = Perc.pos.by.ident(object, features.use = markers, do.plot = FALSE)
    E1 = ExpMat1$ExpMat * ExpMat1$PercMat
    E1 = t(scale(t(E1)))
    # Cluster rows and columns
    test <- heatmap.2(E1, hclustfun = function(x) hclust(x,method="single"))
    if (k==1){
      B1=E1[test$rowInd, rev(test$colInd)]
    } else {
      B1=E1[test$rowInd, test$colInd]
    }
    
    row.max = apply(B1,1,which.max)
    row.ord = names(sort(row.max))
  
   plist1[[k]] = DataHeatmap2(t(B1[row.ord,]), vals.range = c(-1,1), 
                            col.low = "blue", col.mid = "white", col.high = "red", return.plot = TRUE)
   
   # Heatmap for transcriptional clusters
   object = set.all.ident(object, id=id.list[k])
   # Compute Expression matrix
   ExpMat2 = Perc.pos.by.ident(object, features.use = markers, do.plot = FALSE)
   E2 = ExpMat2$ExpMat * ExpMat2$PercMat
   E2 = t(scale(t(E2)))
   # Cluster rows and columns
   B2 = E2[rownames(B1[row.ord,]),]
   col.max = apply(B2,2,which.max)
   col.ord = names(sort(col.max))
   
   plist2[[k]] = DataHeatmap2(t(B2[,col.ord]), vals.range = c(-1,1), 
                              col.low = "blue", col.mid = "white", col.high = "red", return.plot = TRUE)
  
}

dir.create("Fig5")
wid = c(15,18,22,25,30)
hei = c(8,11,14,17,20)
for (k in c(1:5)){
  pdf(paste0("Fig5/",ages[k],"_fatedeterminants.pdf"), w=wid[k],h=hei[k], useDingbats = FALSE)
  plot_grid(plist1[[k]],plist2[[k]], nrow=1)
  dev.off()
  
}

# Overlap between fate markers

fate_markers = list()
for (k in c(1:6)){
  print(k)
  N = N_vec[k]
  eval(parse(text = paste0("object = ", objects[k])))
  if (k < 6){
    object = set.all.ident(object, id="fate_clusters")
  } else {
    object = set.all.ident(object, id="mnn_final")
  }
  print(table(object@ident))
  df = object@de.list$FateClusterDE_MAST
  df$DEratio = df[,3] * df[,7] / ((df[,4]+0.15) * (df[,8]+0.2))
  ExpMat0 = Perc.pos.by.ident(object, features.use = df$gene, do.plot = FALSE)
  E0 = ExpMat0$ExpMat * ExpMat0$PercMat
  
  markers = c()
  for (clust in 1:length(levels(object@ident))){
    df_sub = subset(df, (cluster == clust & diff > 0) & (DEratio > 1.5))
    df_sub = df_sub[order(-df_sub$DEratio),]
    df_sub = subset(df_sub, !(gene %in% markers))
    if (nrow(df_sub) > 0){
      markers = unique(c(markers,df_sub$gene))
      # Find the most specific markers ranked by percentage
    }
    }
    
  markers = markers[!is.na(markers)]
  fate_markers[[ages[k]]] = markers
  
}

# Compute hypergeometric p-values
Pval = matrix(0, nrow=6,ncol=6)
Frac_overlap =  matrix(0, nrow=6,ncol=6)
for (k in c(1:6)){
  for (j in c(2:6)){
    if (k >= j) next
    
    eval(parse(text = paste0("object_k = ", objects[k])))
    eval(parse(text = paste0("object_j = ", objects[j])))
    genes_all = intersect(rownames(object_k@data), rownames(object_j@data))
    genes_k = intersect(fate_markers[[ages[k]]], genes_all)
    genes_j = intersect(fate_markers[[ages[j]]], genes_all)
    m = max(length(genes_j), length(genes_k))
    n = length(genes_all) - m
    k1 = min(length(genes_j), length(genes_k))
    o = length(intersect(genes_k, genes_j))
    Pval[k,j] = -log10(sum(dhyper(o:k1, m, n , k1)))
    Pval[j,k] = Pval[k,j]
    
    Frac_overlap[k,j] = o / length(union(genes_k, genes_j))
    Frac_overlap[j,k] = o/length(union(genes_k, genes_j))
  }
    
  
}

Pval[Pval==Inf] = 300

diag(Pval) = NA
diag(Frac_overlap) = NA

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

colnames(Pval) = ages
colnames(Frac_overlap) = ages
rownames(Pval) = ages
rownames(Frac_overlap) = ages

# Pval
upper_tri <- get_upper_tri(Pval)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat$X1 = factor(melted_cormat$X1, levels = ages)
melted_cormat$X2 = factor(melted_cormat$X2, levels = rev(ages))

# Create a ggheatmap
ggheatmap1 <- ggplot(melted_cormat, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "Red", mid = "orange", 
                       midpoint = 150, limit = c(0,300), space = "Lab", 
                       name="-log10(P)") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap1)

# Pval
upper_tri <- t(get_upper_tri(Frac_overlap))


# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat$X1 = factor(melted_cormat$X1, levels = ages)
melted_cormat$X2 = factor(melted_cormat$X2, levels = rev(ages))

# Create a ggheatmap
ggheatmap2 <- ggplot(melted_cormat, aes(X1, X2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "green", mid = "darkgreen", 
                       midpoint = 0.35, limit = c(0,0.7), space = "Lab", 
                       name="Frac. Overlap") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap2)

pdf("Fig5/OverlapFateDeterminants.pdf",w=14,h=6,useDingbats = FALSE)
plot_grid(ggheatmap1, ggheatmap2)
dev.off()

fate_determinants_all = unique(unlist(fate_markers[1:5]))

genes_all = c()
for (k in c(1:6)){
    eval(parse(text = paste0("object = ", objects[k])))
    genes_all = union(genes_all, rownames(object@data))
}

res_table = list()
for (ont in c("BP","MF","CC")){
  res_table[[ont]] = topGOterms(fg.genes = fate_determinants_all, bg.genes = genes_all, ontology.use = ont)
}

save(list=c("res_table"), file="FateDeterminantsGOanalysis.Rdata")
