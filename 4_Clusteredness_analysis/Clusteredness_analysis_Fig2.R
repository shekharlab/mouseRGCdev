assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

# Functions
source("ClusterednessAnalysisFxns.R")

# load E13
rgcE13 = readRDS("../R_objects/rgcE13_190911_withP56.rds")
rgcE13  = set.all.ident(rgcE13, id="m_merge")
# load E14
rgcE14 = readRDS("../R_objects/rgcE14_190911_withP56.rds")
rgcE14  = set.all.ident(rgcE14, id="m_merge")

# load E16
rgcE16 = readRDS("../R_objects/rgcE16_190911_withP56.rds")
rgcE16  = set.all.ident(rgcE16, id="m_merge")

# load P0
rgcP0 = readRDS("../R_objects/rgcP0_190911_withP56.rds")
rgcP0  = set.all.ident(rgcP0, id="m_merge")

# load P5
rgcP5 = readRDS("../R_objects/rgcP5_190911_withP56.rds")
rgcP5  = set.all.ident(rgcP5, id="mnn_final")

# load atlas
rgc_atlas = readRDS("../R_objects/rgc_atlas_190911_withP56.rds")
rgc_atlas  = set.all.ident(rgc_atlas, id="mnn_final")

objects = c("rgc_atlas","rgcP5", "rgcP0", "rgcE16", "rgcE14","rgcE13")
cluster_ids = c("mnn_final", "m_prune",rep("m_merge", 4))
ages = c("P56","P5","P0","E16","E14","E13")

# load big object
rgc_all = readRDS("../R_objects/rgc_dev_190910.rds")

##########################################
# Within cluster vs cross cluster distance
##########################################

cluster_dists_UMAP = list()
cluster_dists_ligern10 = list()
cluster_dists_ligern20 = list()

for (n in c(1:6)){
  print(ages[n])
  
  # UMAP
  eval(parse(text=paste0("train_object = ", objects[n])))
  train_object = set.all.ident(train_object, id=cluster_ids[n])
  df = as.data.frame(train_object@dr$UMAP@cell.embeddings)
  
  # Standardize each coordinate to -1 and +1
  for (colID in colnames(df)){
    xp = max(df[,colID]); xm = min(df[,colID])
    a = 2/(xp - xm); b = -(xp + xm) / (xp - xm)
    df[,colID] = a*df[,colID] + b
  }
  
  df$cluster = train_object@data.info[,cluster_ids[n]]
  
  res1 = WithinCrossClusterDist(df)
  cluster_dists_UMAP[[ages[n]]]$cross = res1$cross_clust_dist
  cluster_dists_UMAP[[ages[n]]]$self = res1$within_clust_dist
  
}

dist_ratio = rep(0, 6); names(dist_ratio) = rev(ages)
for (n in c(1:6)){
  eval(parse(text=paste0("train_object = ", objects[n])))
  train_object = set.all.ident(train_object, id=cluster_ids[n])
  weights = table(train_object@ident) * length(table(train_object@ident)) / sum(table(train_object@ident))
  x_ratio = weights*cluster_dists_UMAP[[ages[n]]]$cross/cluster_dists_UMAP[[ages[n]]]$self
  dist_ratio[ages[n]] = mean(x_ratio)
}



# Compute validation errors as part of a training paradigm
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
  class_errors_vec[ages[n]] = mean(class_errors[[ages[n]]])
  #class_errors_vec[ages[n]] = sum(class_errors[[ages[n]]] * as.numeric(table(train_object@ident))) / ncol(train_object@data)
  
}

save(list=c("class_errors_vec", "class_errors"), file="Xgboost_class_errors.Rdata")

class_errors_vec = rev(class_errors_vec)


df1 = data.frame(class_errors = class_errors_vec, dist_ratio = 1/dist_ratio)
df1$age = rownames(df1)
pdf("Figs/Clusteredness_metrics_final3.pdf",w=5,h=4, useDingbats = FALSE)
ggplot(df1, aes(x=age, y=class_errors)) + geom_bar(stat="identity", fill="blue") + theme_classic() + xlab("Age") + ylab("Classification Error")
ggplot(df1, aes(x=age, y=dist_ratio)) + geom_bar(stat="identity", fill="red") + theme_classic() + xlab("Age") + ylab("Within Cluster / Cross Cluster")
dev.off()



####### #######
### Compute rao index
################

rgc_all@data.info$all = 1
rgc_all = set.all.ident(rgc_all, id="all")

cluster_ids = c("mnn_final", "mnn_final",rep("m_merge", 4))

rao_index = list()
sd_thresh = c(0.3, 0.5, 0.7, 0.9, 1.1)
l=0
for (iter in c(1:5)){
  print(paste0("iter " = iter))
  var.genes = NB.var.genes(rgc_all,  set.var.genes = FALSE, num.sd = sd_thresh[iter], x.high.cutoff = 30, do.idents = TRUE, x.low.cutoff = 0.0001, rem.mt.rp = FALSE, do.text = FALSE)
  
  rao_index[[iter]] = list()
  rao_index[[iter]]$Ngenes = length(var.genes)
  rao_index[[iter]]$index = c()
  
  for (n in c(1:6)){
    eval(parse(text=paste0("object = ", objects[n])))
    object = set.all.ident(object, id=cluster_ids[n])
    ExpMat = Perc.pos.by.ident(object, features.use = var.genes, do.plot=FALSE)$PercMat
    p_i = table(object@ident) / sum(table(object@ident))
    
    rao_index[[iter]]$index = c(rao_index[[iter]]$index, raoDiversity(t(ExpMat), distfun = "manhattan", p_i = p_i))
    
    if (l == 0){
      df_rao = data.frame(age=ages[n], rao_index = rao_index[[iter]]$index[n], num_genes = rao_index[[iter]]$Ngenes  )
      l=l+1
    } else {
      df_rao_temp = data.frame(age=ages[n], rao_index = rao_index[[iter]]$index[n], num_genes = rao_index[[iter]]$Ngenes  )
      df_rao = rbind(df_rao, df_rao_temp)
    }
  }
}
  
df_rao$age = factor(df_rao$age, levels = rev(ages))
df_rao$rao_index = df_rao$rao_index / df_rao$num_genes

pdf("Figs/Rao_index.pdf", w=6,h=4, useDingbats=FALSE)
ggplot(df_rao, aes(x=age, y=rao_index, group = num_genes)) + geom_point(aes(color=factor(num_genes))) +  geom_line(aes(color=factor(num_genes))) + theme_classic() 
dev.off()


### Shannon and Gini indices

# Diversity
shannonH = rep(0,6); names(shannonH) = rev(levels(rgc_all@data.info$age))
SimpsonI = shannonH
normShannonH = shannonH
normSimpsonI = SimpsonI
for (n in c(1:6)){
  eval(parse(text=paste0("object = ", objects[n])))
  object = set.all.ident(object, id=cluster_ids[n])
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


