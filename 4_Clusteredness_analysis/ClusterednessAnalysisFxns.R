# Packages
library(xgboost)

#' A function to compute the relative tightness of clusters by comparing within cluster diameter to cross-cluster distances
#' 
#' @param coord_ids character vector corresponding to the factors of the reduced dimensional space.
#' @param df data frame with columns corresponding coord_ids and also corresponding to cluster
#' @return a list with elements cross_clust - corresponding to mean cross cluster distances, and within_clust - corresponding to mean within cluster distances
WithinCrossClusterDist = function(df, coord_ids = c("UMAP1","UMAP2")){
  
  if (!("cluster" %in% colnames(df))) stop("df must have a column corresponding to cluster ids")
  
  if (sum(coord_ids %in% colnames(df)) != length(coord_ids) ) stop("df must have columns corresponding to all entries of coord_ids")
  
  if (!is.factor(df$cluster)) df$cluster = factor(df$cluster)

    # Compute centroids
  df_cent = df %>% dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean))
  
  # For each cluster, compute the median within cluster distance and 
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



#' Training an xgboost model for determining classification errors
#'
#' @param train_object0 the scR object used for training
#' @param var.genes character vec, set of features used for training
#' @param do.scale bool, scales the data
#' @param scale.mean numeric vec, user supplied mean vector to use for centering
#' @param scale.var numeric vec, user supplied standard deviation vector for scaling
#'
#' @return list, with elements corresponding to
#' bst_model - the trained xgboost model
#' class_errors - the class specific errors computed from the validation set
#' @export
#'
#' @examples
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




#' Title Computes the rao diversity index
#'
#' @param X matrix with rows as categories and columns as features
#' @param distfun character, distance metric (Default: : `"manhattan"``)
#' @param p_i 
#'
#' @return the value of the Rao diversity index
#' @export
#'
#' @examples
raoDiversity = function(X, distfun = "manhattan", p_i = NULL ){
  
  if (is.null(p_i)){
    p_i = rep(1/nrow(X), nrow(X))
  }
  
  if (length(p_i) != nrow(X)) stop("p_i needs to have the same number of elements as the rows of X")
  if (sum(p_i) != 1) stop("p_i needs to sum to 1")
  
  D = as.matrix(dist(X, method = distfun))
  rao_div = sum (D * (p_i %o% p_i))
  return(rao_div)
}
