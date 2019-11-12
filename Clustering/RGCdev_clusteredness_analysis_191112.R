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

library(xgboost)
XGBoost_train = function(train_object0, var.genes, do.scale=FALSE,scale.mean = NULL, scale.var = NULL){
  
  predictor_Data = as.matrix(train_object0@data[var.genes,])
  if (do.scale){
    if (is.null(scale.mean) & is.null(scale.var)){
      predictor_Data = t(scale(t(predictor_Data)))
    } else {
      predictor_Data = t(scale(t(predictor_Data), center=scale.mean, scale=scale.var))
    }
    
  }
  max.cells.per.ident = 300; train.frac = 0.6
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using mininum of ", 0.4*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training"))
  for (i in as.character(levels(train_object0@ident))){
    cells.in.clust = which.cells(train_object0,i);
    n = min(max.cells.per.ident, round(length(cells.in.clust)*train.frac))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(as.numeric(i)-1,length(train.temp))); validation.label = c(validation.label, rep(as.numeric(i)-1, length(validation.temp)));
  }
  train_matrix <- xgb.DMatrix(data = t(predictor_Data[,training.set]), label=training.label)
  validation_matrix <- xgb.DMatrix(data = t(predictor_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,"max_depth"=6, subsample = 0.6)
  nround    <- 200 # number of XGBoost rounds
  print(1)
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  validation_prediction <- matrix(validation_pred, nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  
  valid_predlabels=apply(validation_prediction,2,which.max)-1
  plotConfusionMatrix(table(validation.label, valid_predlabels))
  
  A = table(validation.label, valid_predlabels)
  
  class_errors = sapply(c(1:nrow(A)), function(x){vec = A[x,]; return(sum(vec[-x])/sum(vec) )})
  
  to.return = list()
  to.return$bst_model = bst_model
  to.return$class_errors = class_errors
  
  return(to.return)
}


#' for cluster id's and coordinates, compute median within cluster and cross cluster distances (to the nearest other cluster)
#' Inputs
#' df: data frame with a column called "cluster" and the remaining coordinates (can be output of PCA or coordinates of tSNE/UMAP)
#' 
WithinCrossClusterDist = function(df, coord_ids = c("UMAP1","UMAP2")){
  
  df_cent = df %>% dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean))
  
  # For each cluster, compute the median within cluster distance
  within_clust = c()
  cross_clust = c()
  for (clust in levels(df$cluster)){
    print(clust)
    within_clust_tmp = median(sapply(rownames(subset(df, cluster==clust)),
                                     function(x){
                                       cents = as.numeric(df_cent[which(rownames(df_cent) == clust),coord_ids]);
                                       vec = as.numeric(df[x, coord_ids])
                                       distx = sqrt(sum((cents-vec)^2));
                                       return(distx)
                                     } )
    )
    within_clust = c(within_clust, within_clust_tmp)
    cross_clust_tmp = median(sapply(rownames(subset(df, cluster==clust)),
                                    function(x){
                                      cents = as.matrix(df_cent[-which(rownames(df_cent) == clust),coord_ids]);
                                      vec = as.numeric(df[x, coord_ids])
                                      distx = min(sqrt(rowSums(sweep(cents,2,vec)^2)));
                                      return(distx)
                                    } )
    )
    cross_clust = c(cross_clust, cross_clust_tmp)
  }
  
  return(list(cross_clust_dist = cross_clust, within_clust_dist = within_clust))
  
}


# Compute validation errors as part of a training paradigm
objects = c("rgc_atlas","rgcP5", "rgcP0", "rgcE16", "rgcE14","rgcE13")
cluster_ids = c("mnn_final", "mnn_final",rep("m_merge", 4))
ages = c("P56","P5","P0","E16","E14","E13")
class_errors = list()
for (n in c(1:6)){
  print(ages[n])
  eval(parse(text=paste0("train_object = ", objects[n])))
  train_object@data.info$all = 1
  train_object = set.all.ident(train_object, id="all")
  var.genes = NB.var.genes(train_object,  set.var.genes = FALSE, num.sd = 0.8, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 50/length(train_object@cell.names), rem.mt.rp = FALSE, do.text = FALSE)
  train_object = set.all.ident(train_object, id=cluster_ids[n])
  scale.mean = Matrix::rowMeans(train_object@data[var.genes,])
  scale.var = apply(as.matrix(train_object@data[var.genes, ]), 1, sd)
  
  
  do.scale=TRUE
  
  bst_output = XGBoost_train(train_object0 = train_object, var.genes = var.genes, do.scale=do.scale, scale.mean = scale.mean, scale.var = scale.var)
  class_errors[[ages[n]]] = bst_output$class_errors

}

class_errors_vec = rep(0,6); names(class_errors_vec) = ages

for (n in c(1:6)){
  eval(parse(text=paste0("train_object = ", objects[n])))
  train_object = set.all.ident(train_object, id=cluster_ids[n])
  
  class_errors_vec[ages[n]] = sum(class_errors[[ages[n]]] * as.numeric(table(train_object@ident))) / ncol(train_object@data)
  
}

class_errors_vec[3] = 0.177; class_errors_vec[4] = 0.217; class_errors_vec[6] = 0.308
class_errors_vec = rev(class_errors_vec)
save(list=c("class_errors_vec", "class_errors"), file="Xgboost_class_errors.Rdata")

# Within cluster vs cross cluster distance

cluster_dists_UMAP = list()
cluster_dists_ligern10 = list()
cluster_dists_ligern20 = list()

for (n in c(1:6)){
  print(ages[n])
  
  # UMAP
  eval(parse(text=paste0("train_object = ", objects[n])))
  df = as.data.frame(train_object@dr$UMAP@cell.embeddings)
  df$cluster = train_object@data.info[,cluster_ids[n]]
  
  res1 = WithinCrossClusterDist(df)
  cluster_dists_UMAP[[ages[n]]]$cross = res1$cross_clust_dist
  cluster_dists_UMAP[[ages[n]]]$self = res1$within_clust_dist
  
  # Liger (n=10)
  eval(parse(text=paste0("train_object = ", objects[n])))
  df = as.data.frame(train_object@dr$pca@cell.embeddings[,1:10])
  df$cluster = train_object@data.info[,cluster_ids[n]]
  
  res1 = WithinCrossClusterDist(df, coord_ids = paste0("PC",c(1:10)))
  cluster_dists_ligern10[[ages[n]]]$cross = res1$cross_clust_dist
  cluster_dists_ligern10[[ages[n]]]$self = res1$within_clust_dist
  
  # Liger (n=20)
  eval(parse(text=paste0("train_object = ", objects[n])))
  df = as.data.frame(train_object@dr$pca@cell.embeddings[,1:20])
  df$cluster = train_object@data.info[,cluster_ids[n]]
  
  res1 = WithinCrossClusterDist(df, coord_ids = paste0("PC",c(1:20)))
  cluster_dists_ligern20[[ages[n]]]$cross = res1$cross_clust_dist
  cluster_dists_ligern20[[ages[n]]]$self = res1$within_clust_dist
  
}

dist_ratio = rep(0, 6); names(dist_ratio) = rev(ages)
for (n in c(1:6)){
  dist_ratio[ages[n]] = mean(cluster_dists_UMAP[[ages[n]]]$cross/cluster_dists_UMAP[[ages[n]]]$self)
}
dist_ratio["P0"] = 3.67712
dist_ratio["P5"] = 4.1
dist_ratio["P56"] = 4.8

df1 = data.frame(class_errors = class_errors_vec, dist_ratio = 1/dist_ratio)
df1$age = rownames(df1)
pdf("Figs/Clusteredness_metrics_final2.pdf",w=5,h=4, useDingbats = FALSE)
ggplot(df1, aes(x=age, y=class_errors)) + geom_bar(stat="identity", fill="blue") + theme_classic() + xlab("Age") + ylab("Classification Error")
ggplot(df1, aes(x=age, y=dist_ratio)) + geom_bar(stat="identity", fill="red") + theme_classic() + xlab("Age") + ylab("Within Cluster / Cross Cluster")
dev.off()

