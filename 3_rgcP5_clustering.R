assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

rgcP5 = readRDS("../R_objects/rgcP5_190729.rds")
rgcP5 = set.all.ident(rgcP5, id="orig")
DimPlot(rgcP5, reduction.use = "UMAP")
rgcP5 = set.all.ident(rgcP5, id="m")
DimPlot(rgcP5, reduction.use = "UMAP")

cells.use=colnames(rgcP5@data)
batchname = sapply(cells.use,getStat1)
batchid = rep("Batch5",length(cells.use))
batchid[grep("S1$",batchname)] = "Batch1"
batchid[grep("S2$",batchname)] = "Batch1"
batchid[grep("S3$",batchname)] = "Batch1"
batchid[grep("S4$",batchname)] = "Batch2"
batchid[grep("S5$",batchname)] = "Batch2"
batchid[grep("S6$",batchname)] = "Batch2"
batchid[grep("S7$",batchname)] = "Batch2"
batchid[grep("S8$",batchname)] = "Batch2"
batchid[grep("S13$",batchname)] = "Batch2"
batchid[grep("S9$",batchname)] = "Batch3"
batchid[grep("S10$",batchname)] = "Batch3"
batchid[grep("S11$",batchname)] = "Batch3"
batchid[grep("S12$",batchname)] = "Batch3"
batchid[grep("S14$",batchname)] = "Batch4"
batchid[grep("S15$",batchname)] = "Batch4"
batchid[grep("S16$",batchname)] = "Batch4"
rgcP5@data.info$batch = factor(batchid)

# Ligerize
rgcP5@data.info$all = 1
rgcP5 = set.all.ident(rgcP5, id="all")
var.genes = NB.var.genes(rgcP5,  set.var.genes = FALSE, num.sd = 1, x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0005, rem.mt.rp = FALSE, do.text = FALSE)
rgcP5@var.genes = var.genes

rgcP5 = Ligerize(rgcP5, batch_id = "batch", var.genes = var.genes, do.umap = TRUE)
saveRDS(rgcP5, file="../R_objects/rgcP5_190729.rds")

rgcP5 = set.all.ident(rgcP5, id="orig")
DimPlot(rgcP5, reduction.use="UMAP")

rgcP5 = set.all.ident(rgcP5, id="m")
DimPlot(rgcP5, reduction.use="UMAP", do.label=TRUE)

p=Perc.pos.by.ident(rgcP5,c("Gap43","Rtn1","Sncg","Rbpms","Nefl","Pou4f1","Thy1", "Slc17a6","Onecut2", "Pax6", "Tfap2b","Tfap2a","Gad1","Slc6a9", "C1qa","H19","Gngt2", "Ccnd1","Hmgb2","Fgf15","Hes5","Id3","Btg1","Snap25"), max.val.exp = 5, alpha.use=0.7, return.plot = TRUE)
p

a=find.markers(rgcP5, 28, test.use="MAST", max.cells.per.ident = 1000)


rgcP5 = prune.clust(rgcP5, remove.clust = c(28,31,32,33))


ExpMat = Avg.by.ident(rgcP5, features.use = rgcP5@var.genes,return.scaled=FALSE)
a = apply(ExpMat,1, function(x) if(max(x) < 0.3){return(NA)} else {return(mean(sort(x,decreasing=TRUE)[1:7]) / mean(sort(x,decreasing=FALSE)[1:7]))})
a=a[!is.na(a)]
var.genes = names(a)[a>5]
library(lsa)
rgcP5 = buildClusterTree(rgcP5, genes.use=rgcP5@var.genes, regress.use = FALSE,linkage.method = "average",dist.fun = "cosine")

levels(rgcP5@ident) = 1:length(levels(rgcP5@ident))
rgcP5@data.info$m_merge = rgcP5@ident
rgcP5 = find.all.markers(rgcP5, test.use="MAST", max.cells.per.ident = 2000)
saveRDS(rgcP5, file="../R_objects/rgcP5_190729.rds")

# Read old p5RGC object
rgcP5_old = readRDS("../R_objects/p5RGC_cluster_final_2018-05-16.rds")
rgcP5_old = rgcP5_old$p5RGC_clean
old_ident = rgcP5_old@ident
names(old_ident) = gsub("p5RGC","P5Cd90RGC",names(old_ident))
names(old_ident) = gsub("10x","10X",names(old_ident))
cells.common = intersect(names(old_ident), rgcP5@cell.names)

train_object = subsetData(rgcP5, cells.use = cells.common)
train_object@data.info$mnn_final = old_ident[train_object@cell.names]
train_object = set.all.ident(train_object, id="mnn_final")

# Transfer labels by first training a classifier
bst = XGBoost_train(train_object0 = train_object, var.genes = train_object@var.genes, do.scale=TRUE)
numberOfClasses <- 38
test_Data = rgcP5@data[rgcP5@var.genes,]
test_Data = t(scale(t(test_Data)))
test_matrix <- xgb.DMatrix(data = t(test_Data))
test_pred <- predict(bst$bst_model, newdata = test_matrix)
test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                          ncol=length(test_pred)/numberOfClasses)
test_pred_margins = apply(test_prediction,2,max)
test_predlabels=apply(test_prediction,2,which.max)
names(test_predlabels) = colnames(rgcP5@data)

rgcP5@data.info$mnn_final = factor(test_predlabels[rgcP5@cell.names])
rgcP5 = set.all.ident(rgcP5, id="mnn_final")
saveRDS(rgcP5, file="../Objects/rgcP5_190729.rds")

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
  to.return$Conf = A
  to.return$true_label = validation.label
  to.return$assigned_label = valid_predlabels
  
  return(to.return)
}
