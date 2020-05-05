assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

source("../scripts/AllFxns.R")

# load E13
rgcE13 = readRDS("../R_objects/rgcE13_190729.rds")
# load E14
rgcE14 = readRDS("../R_objects/rgcE14_190729.rds")
# load E16
rgcE16 = readRDS("../R_objects/rgcE16_190729.rds")
# load P0
rgcP0 = readRDS("../R_objects/rgcP0_190729.rds")
# load P5
rgcP5 = readRDS("../R_objects/rgcP5_190729.rds")

# load atlas
rgc_atlas = readRDS("../R_objects/rgc_atlas_190729.rds")

# load RGC dev
rgc_dev = readRDS("../R_objects/rgc_dev_190729.rds")

var.genes = rgc_atlas@var.genes

# Find the common variable genes
objects = c("rgcP5", "rgcP0", "rgcE16", "rgcE14","rgcE13")
cluster_ids = c("mnn_final",rep("m_merge", 4))
for (n in c(1:5)){
  eval(parse(text=paste0("object = ", objects[n])))
var.genes = var.genes[var.genes %in% rownames(object@data)]
}

# Train model on atlas
rgc_dev = set.all.ident(rgc_dev, id="orig")
rgc_atlasx = subsetData(rgc_dev, cells.use = which.cells(rgc_dev, grep("^aRGC", levels(rgc_dev@ident), value=TRUE)))
rgc_atlasx@data.info$m = rgc_atlas@data.info[rownames(rgc_atlasx@data.info), "mnn_final"]
rgc_atlasx = set.all.ident(rgc_atlasx, id="m")
scale.mean = Matrix::rowMeans(rgc_atlasx@data[var.genes,])
scale.var = apply(as.matrix(rgc_atlasx@data[var.genes, ]), 1, sd)
do.scale=TRUE
xgboost_atlas = XGBoost_train(train_object0 = rgc_atlasx, var.genes = var.genes, do.scale=do.scale, scale.mean = scale.mean, scale.var = scale.var, min.val = -5, max.val = 7)
imp_atlas = as.data.frame(xgb.importance(var.genes, model = xgboost_atlas$bst_model))

save(list=c("xgboost_atlas", "var.genes","imp_atlas"), file="Xgboost_for_classifying_DevRGCs_02182020.Rdata")
load("Xgboost_for_classifying_DevRGCs_02182020.Rdata")


######################################################################
####### First Classify Each age to atlas ############################
###################################################################

atlas_xgboost_res = list()
cluster_ids = c("mnn_final",rep("m_merge", 4))
ages = c("P5","P0","E16","E14","E13")
atlas_max = apply(rgc_atlasx@data[var.genes,],1,max)
xgboost_dev_classification = list()

for (n in c(1:5)){
  
  xgboost_dev_classification[[ages[n]]] = list()
  test_object = subsetData(rgc_dev, cells.use = grep(paste0("^",ages[n]), rgc_dev@cell.names, value=TRUE))
  eval(parse(text = paste0("test_object@data.info[,'",cluster_ids[n],"'] = ", objects[n], "@data.info[test_object@cell.names, '", cluster_ids[n],"']" )))
  test_object = set.all.ident(test_object, id= cluster_ids[n])
  do.scale=TRUE
  xgboost_test = XGBoost_train(train_object0 = test_object, var.genes = var.genes, do.scale=do.scale, scale.mean = scale.mean, scale.var = scale.var,  min.val = -5, max.val = 7)
  imp_test = as.data.frame(xgb.importance(var.genes, model = xgboost_test$bst_model))

  # Intersection of the top 500 genes
  genes.use = intersect(imp_atlas$Feature[1:500], imp_test$Feature[1:500])
  
  # Retrain atlas
  do.scale=TRUE
  xgboost_atlas_new = XGBoost_train(train_object0 = rgc_atlasx, var.genes = genes.use, do.scale=do.scale, scale.mean = scale.mean[genes.use], scale.var = scale.var[genes.use],  min.val = -5, max.val = 7)
  test_Data = as.matrix(test_object@data[genes.use,])
  #if (do.scale) test_Data = t(scale(t(test_Data), scale =  scale.var[genes.use], center=scale.mean[genes.use]))
  if (do.scale) test_Data = t(scale(t(test_Data)))
  test_Data = minmax(test_Data, min = -5, max = 7)
  
  numberOfClasses <- 45
  test_matrix <- xgb.DMatrix(data = t(test_Data))
  test_pred <- predict(xgboost_atlas_new$bst_model, newdata = test_matrix)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses)
  
  test_pred_margins = apply(test_prediction,2,max)
  test_predlabels=apply(test_prediction,2,which.max)
  names(test_predlabels) = colnames(test_object@data)
  test_pred_names = levels(rgc_atlas@ident)[test_predlabels]
  names(test_pred_names) = names(test_predlabels)
  
  xgboost_dev_classification[[ages[n]]]$atlas_labels = test_predlabels
  xgboost_dev_classification[[ages[n]]]$orig_cluster = test_object@ident
}


save(list=c("xgboost_atlas", "var.genes","imp_atlas", "xgboost_dev_classification"), file="Xgboost_for_classifying_DevRGCs_02182020.Rdata")

load("Xgboost_classifcation_each_stage_02182020.Rdata")

######################################################################
####### Next, Classify Each age to next time point ############################
###################################################################

objects = c("rgcE13","rgcE14","rgcE16","rgcP0","rgcP5","rgc_atlas")
cluster_ids = c(rep("m_merge", 4),"mnn_final","mnn_final")
ages = rev(c("aRGC", "P5","P0","E16","E14","E13"))
xgboost_dev_classification_next = list()
for (n in c(1:6)){
  print(n)
  if (n < 6){
    test_object = subsetData(rgc_dev, cells.use = grep(paste0("^",ages[n]), rgc_dev@cell.names, value=TRUE))
    train_object = subsetData(rgc_dev, cells.use = grep(paste0("^",ages[n+1]), rgc_dev@cell.names, value=TRUE))
    eval(parse(text = paste0("train_object@data.info[,'",cluster_ids[n+1],"'] = ", objects[n+1], "@data.info[train_object@cell.names, '", cluster_ids[n+1],"']" )))
    eval(parse(text = paste0("test_object@data.info[,'",cluster_ids[n],"'] = ", objects[n], "@data.info[test_object@cell.names, '", cluster_ids[n],"']" )))
    
    train_object = set.all.ident(train_object, id= cluster_ids[n+1])
    
  } else {
    test_object = subsetData(rgc_dev, cells.use = grep(paste0("^",ages[n]), rgc_dev@cell.names, value=TRUE))
    train_object = subsetData(rgc_dev, cells.use = grep(paste0("^",ages[n-1]), rgc_dev@cell.names, value=TRUE))
    
    eval(parse(text = paste0("train_object@data.info[,'",cluster_ids[n-1],"'] = ", objects[n-1], "@data.info[train_object@cell.names, '", cluster_ids[n-1],"']" )))
    eval(parse(text = paste0("test_object@data.info[,'",cluster_ids[n],"'] = ", objects[n], "@data.info[test_object@cell.names, '", cluster_ids[n],"']" )))
    
    train_object = set.all.ident(train_object, id= cluster_ids[n-1])
  }
  
  
  xgboost_dev_classification_next[[ages[n]]] = list()
  
  # Train training object 
  scale.mean = Matrix::rowMeans(train_object@data[var.genes,])
  scale.var = apply(as.matrix(train_object@data[var.genes, ]), 1, sd)
  xgboost_train = XGBoost_train(train_object0 = train_object, var.genes = var.genes, do.scale=TRUE, scale.mean = scale.mean, scale.var = scale.var,   min.val = -5, max.val = 7)
  imp_train = as.data.frame(xgb.importance(var.genes, model = xgboost_train$bst_model))
  
  # Train test object
  test_object = set.all.ident(test_object, id= cluster_ids[n])
  xgboost_test = XGBoost_train(train_object0 = test_object, var.genes = var.genes, do.scale=TRUE, scale.mean = scale.mean, scale.var = scale.var,   min.val = -5, max.val = 7)
  imp_test = as.data.frame(xgb.importance(var.genes, model = xgboost_test$bst_model))
  
  # Intersection of the top 500 genes
  genes.use = intersect(imp_train$Feature[1:500], imp_test$Feature[1:500])
  genes.use = genes.use[!is.na(genes.use)]
  # Retrain training object
  do.scale=TRUE
  xgboost_train_new = XGBoost_train(train_object0 = train_object, var.genes = genes.use, do.scale=do.scale, scale.mean = scale.mean[genes.use], scale.var = scale.var[genes.use],   min.val = -5, max.val = 7)
  test_Data = as.matrix(test_object@data[genes.use,])
  if (do.scale) test_Data = t(scale(t(test_Data)))
  test_Data = minmax(test_Data, min = -5, max = 7)
  
  numberOfClasses <- length(levels(train_object@ident))
  test_matrix <- xgb.DMatrix(data = t(test_Data))
  test_pred <- predict(xgboost_train_new$bst_model, newdata = test_matrix)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses)
  
  test_pred_margins = apply(test_prediction,2,max)
  test_predlabels=apply(test_prediction,2,which.max)
  names(test_predlabels) = colnames(test_object@data)
  test_pred_names = levels(rgc_atlas@ident)[test_predlabels]
  names(test_pred_names) = names(test_predlabels)
  
  xgboost_dev_classification_next[[ages[n]]]$new_cluster = test_predlabels
  xgboost_dev_classification_next[[ages[n]]]$orig_cluster = test_object@ident
}

save(list=c("ari_xgboost","ari_xgboost_next", "xgboost_dev_classification", "xgboost_dev_classification_next"),
     file="Xgboost_classifcation_each_stage_02182020.Rdata")

load("Xgboost_classifcation_each_stage_02182020.Rdata")

#########################################################
###### Compute ARI and NCE scores ######################
#########################################################


ari_xgboost = rep(0,5); 
ages = rev(c("P5","P0","E16","E14","E13"))
names(ari_xgboost) = ages
ari_xgboost_rand = ari_xgboost
library(mclust)
for (n in c(1:5)){
  ari_xgboost[ages[n]] = adjustedRandIndex(xgboost_dev_classification[[ages[n]]]$atlas_labels, xgboost_dev_classification[[ages[n]]]$orig_cluster)
  ari_xgboost_rand[ages[n]] = abs(adjustedRandIndex(sample(xgboost_dev_classification[[ages[n]]]$atlas_labels), xgboost_dev_classification[[ages[n]]]$orig_cluster))
}

# Compute validation errors as part of a training paradigm


ari_xgboost_next = rep(0,5); 
ages = rev(c("P56", "P5","P0","E16","E14","E13"))
names(ari_xgboost_next) = ages[1:5]
ari_xgboost_next_rand = ari_xgboost_next

library(mclust)
for (n in c(1:5)){
  ari_xgboost_next[ages[n]] = adjustedRandIndex(xgboost_dev_classification_next[[ages[n]]]$new_cluster, xgboost_dev_classification_next[[ages[n]]]$orig_cluster)
  ari_xgboost_next_rand[ages[n]] = abs(adjustedRandIndex(sample(xgboost_dev_classification_next[[ages[n]]]$new_cluster), xgboost_dev_classification_next[[ages[n]]]$orig_cluster))

}


df = data.frame(ARI_atlas = rev(ari_xgboost) ,ARI_next = rev(ari_xgboost_next),  age = ages[1:5] )
df$age = rownames(df)
df_melt = melt(df)
colnames(df_melt)[c(2,3)] = c("Method","ARI")
dfARI = df

dir.create("Figs")
pdf("Figs/ARI_Xgboost_dev_final.pdf",w=6,h=4)
ggplot(df_melt, aes(x=age, y=ARI)) + geom_bar(stat="identity", aes(fill=factor(Method)),color="black", position="dodge") + theme_classic() + xlab("Age") + ylab("Adjusted Rand Index") + scale_fill_manual(values = c("cyan","lightslateblue"))
dev.off()


cluster_ids = rev(c("mnn_final",rep("m_merge", 4)))
ages = rev(c("P56", "P5","P0","E16","E14","E13"))
cond_entropy_atlas = rep(0,5); names(cond_entropy_atlas) = ages[1:5]
cond_entropy_next = cond_entropy_atlas
library(infotheo)
for (n in c(1:5)){
  e1 = condentropy(xgboost_dev_classification[[ages[n]]]$orig_cluster,xgboost_dev_classification[[ages[n]]]$atlas_labels)
  e1 = e1 / entropy(xgboost_dev_classification[[ages[n]]]$orig_cluster)
  e2 = condentropy(xgboost_dev_classification_next[[ages[n]]]$orig_cluster,xgboost_dev_classification_next[[ages[n]]]$new_cluster)
  e2 = e2 / entropy(xgboost_dev_classification_next[[ages[n]]]$orig_cluster)
  
  #e1 = condentropy(xgboost_dev_classification[[ages[n]]]$atlas_labels,xgboost_dev_classification[[ages[n]]]$orig_cluster)
  #e1 = e1 / entropy(xgboost_dev_classification[[ages[n]]]$atlas_labels)
  #e2 = condentropy(xgboost_dev_classification_next[[ages[n]]]$new_cluster,xgboost_dev_classification_next[[ages[n]]]$orig_cluster)
  #e2 = e2 / entropy(xgboost_dev_classification_next[[ages[n]]]$new_cluster)
  
  cond_entropy_atlas[ages[n]] = e1
  cond_entropy_next[ages[n]] = e2
}


df = data.frame(NormCE_atlas = cond_entropy_atlas , NormCE_next = cond_entropy_next,  age = ages[1:5] )
dfNCE = df
df_melt = melt(df)
colnames(df_melt)[c(2,3)] = c("Method","NormCE")
pdf("Figs/NormCondEntropy_Xgboost_dev_final.pdf",w=6,h=4)
ggplot(df_melt, aes(x=age, y=NormCE)) + geom_bar(stat="identity", aes(fill=factor(Method)),color="black", position="dodge") + theme_classic() + xlab("Age") + ylab("Normalized Conditional Entropy (Xgboost)") + scale_fill_manual(values = c("cyan","lightslateblue"))
dev.off()

# Confusion matrices

A=table(xgboost_dev_classification_next[["E13"]]$new_cluster, xgboost_dev_classification_next[["E13"]]$orig_cluster)
test <- heatmap.2(A,hclustfun = function(x) hclust(x,method="single"))
B=A[test$rowInd, test$colInd]
row.max = apply(B,1,which.max)
row.ord = names(sort(row.max))
pdf("ConfusionMatrixNext_E13_E14.pdf",w=6,h=5,useDingbats = FALSE)
plotConfusionMatrix(B[row.ord,], col.low = "white", col.high = "darkblue", ylab.use = "P56 labels", xlab.use = "P5 clusters")
dev.off()

# P5 vs P56
A=table(xgboost_dev_classification_next[["P5"]]$new_cluster, xgboost_dev_classification_next[["P5"]]$orig_cluster)
test <- heatmap.2(A, hclustfun = function(x) hclust(x,method="single"))
B=A[test$rowInd, test$colInd]
row.max = apply(B,1,which.max)
row.ord = names(sort(row.max))
pdf("ConfusionMatrixNext_P5_P56.pdf",w=10,h=8,useDingbats = FALSE)
plotConfusionMatrix(B[row.ord,], col.low = "white", col.high = "darkblue", ylab.use = "P56 labels", xlab.use = "P5 clusters")
dev.off()

# Plot all confusion matrices for the supplement
ages = rev(c("P56", "P5","P0","E16","E14","E13"))
max_perc_vals = c(25, 40, 60,80,100)
plist = list()
cluster_groups = list()

for (n in c(1:5)){
  A=table(xgboost_dev_classification[[ages[n]]]$atlas_labels, xgboost_dev_classification[[ages[n]]]$orig_cluster)
  A_scale = scale(A, center=FALSE, scale=colSums(A))
  test <- heatmap.2(A_scale, hclustfun = function(x) hclust(x,method="single"))
  B=A[test$rowInd, test$colInd]
  row.max = apply(B,1,which.max)
  row.ord = names(sort(row.max))
  ARI = dfARI[ages[n],1]
  NCE = dfNCE[ages[n],1]
  plist[[n]]=plotConfusionMatrix(B[row.ord,], row.scale=TRUE,col.low = "white", col.high = "darkblue", ylab.use = "P56 labels", xlab.use = paste0(ages[n]," clusters"), plot.return = TRUE, max.perc =100) + ggtitle(paste0("NCE = ", round(NCE,2),", ARI =", round(ARI,2) ))
  
  # With next
  A=table(xgboost_dev_classification_next[[ages[n]]]$new_cluster, xgboost_dev_classification_next[[ages[n]]]$orig_cluster)
  A_scale = scale(A, center=FALSE, scale=colSums(A))
  test <- heatmap.2(A_scale, hclustfun = function(x) hclust(x,method="single"))
  B=A[test$rowInd, test$colInd]
  row.max = apply(B,1,which.max)
  row.ord = names(sort(row.max))
  ARI = dfARI[ages[n],2]
  NCE = dfNCE[ages[n],2]
  plist[[n+5]]=plotConfusionMatrix(B[row.ord,], row.scale=TRUE, col.low = "white", col.high = "darkblue", ylab.use = paste0(ages[n+1]," labels"), xlab.use = paste0(ages[n]," clusters"), plot.return = TRUE, max.perc = 100) + ggtitle(paste0("NCE = ", round(NCE,2),", ARI =", round(ARI,2) ))
  
  A_scale = A_scale[test$rowInd, test$colInd][row.ord,]
  # Define cluster groups
  cluster_groups[[ages[n]]] = FindClusterGroups(X=A_scale, corr_thresh = 0.6)
  
}

pdf("Figs/ConfusionMatrices_clusters_Atlas_final_rowscale.pdf",w=35,h=8, useDingbats = FALSE)
plot_grid(plotlist = plist[1:5], nrow=1, rel_widths = c(1,1.2,1.4,1.6,2))
dev.off()

pdf("Figs/ConfusionMatrices_clusters_Next_final_rowscale.pdf",w=35,h=8, useDingbats = FALSE)
plot_grid(plotlist = plist[6:10], nrow=1, rel_widths = c(1,1.2,1.4,1.6,2), rel_heights = c(1,1.2,1.4,1.6,2))
dev.off()


### Collapse into groups

# Define cell groups
objects = c("rgcE13","rgcE14","rgcE16","rgcP0","rgcP5","rgc_atlas")
cluster_ids = c(rep("m_merge", 4),"mnn_final","mnn_final")
ages = rev(c("P56", "P5","P0","E16","E14","E13"))

plist = list()

cluster_groups[[5]][26:34] = cluster_groups[[5]][25:33]
cluster_groups[[5]][24] = list(10)
cluster_groups[[5]][25] = list(29)

# Xgboost vs. OT confusion matrices
for (n in c(1:5)){
  
  eval(parse(text=paste0("object_now = ", objects[n])))
  eval(parse(text=paste0("object_next = ", objects[n+1])))
  
  object_now = set.all.ident(object_now, id=cluster_ids[n])
  object_next = set.all.ident(object_next, id=cluster_ids[n])
  
  # Xgboost 
  A=table(xgboost_dev_classification_next[[ages[n]]]$new_cluster, xgboost_dev_classification_next[[ages[n]]]$orig_cluster)
  group_vec = getGroupOrder(cluster_groups[[ages[n]]])
  A = t(rowsum(t(A), group=group_vec))
  colnames(A) = getGroupNames(cluster_groups[[ages[n]]])
  
  if (n < 5){
  group_vec = getGroupOrder(cluster_groups[[ages[n+1]]])
  A = rowsum(A, group=group_vec)
  rownames(A) = getGroupNames(cluster_groups[[ages[n+1]]])
} 
  
  A_scale = scale(A, center=FALSE, scale=colSums(A))
  row.max = apply(A_scale,1,which.max)
  row.ord = names(sort(row.max))
  
  plist[[n]]=plotConfusionMatrix(A[row.ord,], row.scale=TRUE, col.low = "white", col.high = "darkblue", ylab.use = paste0(ages[n+1]," groups"), xlab.use = paste0(ages[n]," groups"), plot.return = TRUE, max.perc = 0.9*max(A_scale)*100, x.lab.rot = TRUE, alpha.use = 0.8) 
  
  # OT
  
  M = read.table(paste0("../../OT_analysis/OT_uncorr_vargenes_eps005/tmaps/wot_eps005_",ages[n+1],"_traj.txt_trajectory.txt"), 
                 header=TRUE, row.names = 1, sep="\t")
  M = M[object_now@cell.names,]
  
  # Sum by cluster
  B = rowsum(M, group=object_now@data.info[,cluster_ids[n]])
  
  
  # Collapse into groups
  B = rowsum(B, getGroupOrder(cluster_groups[[ages[n]]]))
  rownames(B) = getGroupNames(cluster_groups[[ages[n]]])
  B = t(B)
  
  if (n < 5){
    B = rowsum(B, getGroupOrder(cluster_groups[[ages[n+1]]]))
    rownames(B) = getGroupNames(cluster_groups[[ages[n+1]]])
  } else {
    rownames(B) = c(1:45) 
  }
  
  
  plist[[n+5]] = plotConfusionMatrix(B[row.ord,], row.scale=TRUE, col.low = "white", col.high = "darkred", ylab.use = paste0(ages[n+1]," groups"), xlab.use = paste0(ages[n]," groups"), plot.return = TRUE, max.perc = 0.9*max(A_scale)*100, x.lab.rot = TRUE, alpha.use = 0.8) 
  
  
}

save(list=c("cluster_groups"), file="Cluster_groups.Rdata")

pdf("Figs/E13_E14_cluster_groups.pdf",w=10,h=3, useDingbats = FALSE)
plot_grid(plist[[1]], plist[[6]])
dev.off()

pdf("Figs/E14_E16_cluster_groups.pdf",w=10,h=5, useDingbats = FALSE)
plot_grid(plist[[2]], plist[[7]])
dev.off()

pdf("Figs/E16_P0_cluster_groups.pdf",w=12,h=6, useDingbats = FALSE)
plot_grid(plist[[3]], plist[[8]])
dev.off()

pdf("Figs/P0_P5_cluster_groups.pdf",w=14,h=7, useDingbats = FALSE)
plot_grid(plist[[4]], plist[[9]])
dev.off()


pdf("Figs/P5_P56_cluster_groups.pdf",w=16,h=8, useDingbats = FALSE)
plot_grid(plist[[5]], plist[[10]])
dev.off()



# Cartoon
Sp = zeros(8,4)
Sp[c(1,2),1] = c(1,1)
Sp[c(3,4,5),2] = c(1,1,1)
Sp[c(6),3] = c(1)
Sp[c(7,8),4] = c(1,1)

pdf("Figs/Cartoon_specified.pdf",w=4,h=4, useDingbats = FALSE)
plotConfusionMatrix(Sp, row.scale=TRUE,col.low = "white", col.high = "darkblue", max.perc=100)
dev.off()

Sp = matrix(runif(32), nrow=8)

pdf("Figs/Cartoon_nonspecified.pdf",w=4,h=4, useDingbats = FALSE)
plotConfusionMatrix(Sp, row.scale=TRUE,col.low = "white", col.high = "darkblue", max.perc=100)
dev.off()


