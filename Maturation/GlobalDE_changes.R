assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")


rgc_dev = readRDS("../Objects/rgc_dev_190910.rds")
rgc_dev = set.all.ident(rgc_dev, id="age")
cells.use = subsample_by_ident(rgc_dev@ident, downsample.frac = 1, max.cells.per.ident = 4000)
rgc_dev = subsetData(rgc_dev, cells.use = cells.use)

##################################################
########### Find global genes that change significantly
##################################################
rgc_dev@data.info$all = 1
rgc_dev = set.all.ident(rgc_dev, id = "all")
var.genes = NB.var.genes(rgc_dev,  set.var.genes = FALSE, num.sd = 0.3, x.high.cutoff = 50, do.idents = TRUE, x.low.cutoff = 0.0005, rem.mt.rp = FALSE, do.text = FALSE)
var.genes = unique(c(var.genes, rgc_dev@var.genes))
rgc_dev = set.all.ident(rgc_dev, id="age")
ExpMat_global = Perc.pos.by.ident(rgc_dev, features.use = var.genes, do.plot = FALSE)
orig_ident = rgc_dev@ident


# Randomize
rgc_dev@ident = sample(orig_ident); names(rgc_dev@ident) = names(orig_ident)
ExpMat_rand_global = Perc.pos.by.ident(rgc_dev, features.use = var.genes, do.plot = FALSE)
rgc_dev@ident = orig_ident

size_fac = (rgc_dev@data.info %>% dplyr::group_by(age) %>% dplyr::summarise_all(mean))$nTranscripts
size_fac = size_fac / mean(size_fac)
change_fac = sapply(var.genes, function(j){
  if (max(ExpMat_global$PercMat[j,]) < 0.2) { return(0) }
  x = ExpMat_global$PercMat[j,] * ExpMat_global$ExpMat[j,] * size_fac
  if ((max(x) - min(x)) / min(x) < 0.3){ return(0) }
  
  a=ExpMat_global$PercMat[j,] * ExpMat_global$ExpMat[j,]; b=ExpMat_rand_global$PercMat[j,] * ExpMat_rand_global$ExpMat[j,];
  #return(sum((a-b)^2) / sqrt(mean(a^2)*mean(b^2)))
  return(sum((a-b)^2) / sqrt(mean((a+0.1)^2)*mean((b+0.1)^2)))
  #return(sum((a-b)^2) * sqrt(sum(diff(a,1)^2))/sqrt(sum(diff(b,1)^2)))
})

change_fac = sort(change_fac, decreasing=TRUE)

# Randomized change fac
change_fac_rand = sapply(var.genes, function(j){
  if (max(ExpMat_global$PercMat[j,]) < 0.2) { return(0) }
  x = ExpMat_global$PercMat[j,] * ExpMat_global$ExpMat[j,] * size_fac
  if ((max(x) - min(x)) / min(x) < 0.3){ return(0) }
  
  a=ExpMat_rand2_global$PercMat[j,] * ExpMat_rand2_global$ExpMat[j,]; b=ExpMat_rand_global$PercMat[j,] * ExpMat_rand_global$ExpMat[j,];
  #return(sum((a-b)^2) / sqrt(mean(a^2)*mean(b^2)))
  return(sum((a-b)^2) / sqrt(mean((a+0.1)^2)*mean((b+0.1)^2)))
  #return(sum((a-b)^2) * sqrt(sum(diff(a,1)^2))/sqrt(sum(diff(b,1)^2)))
})

change_fac_rand = sort(change_fac_rand, decreasing=TRUE)

time_max = sapply(var.genes, function(j){
  return(colnames(ExpMat_global$PercMat)[which.max(ExpMat_global$PercMat[j,] * ExpMat_global$ExpMat[j,] * size_fac)])
})
time_max = time_max[names(change_fac)]

genes.use = names(change_fac)[change_fac > 0.7]
A = ExpMat_global$ExpMat[genes.use,] * ExpMat_global$PercMat[genes.use,] * t(replicate(length(genes.use), size_fac))
A = t(scale(t(A))) 
wss = sapply(1:25,function(k){ kmeans(A,k, nstart=10, iter.max = 50)$tot.withinss})
secondDerivative = c()
for (k in 2:24){
  secondDerivative = c(secondDerivative, wss[k+1] + wss[k-1] - 2*wss[k])
}

nclust = which.min(secondDerivative) + 1
nclust = 6
set.seed(1234)
a=kmeans(A, nclust, nstart = 20)

# Find cluster order
cluster_order = rep(0, nclust)
names(cluster_order) = c(1:nclust)
for (k in 1:nclust){
  genes.use = names(a$cluster)[a$cluster == k]
  meanA = apply(A[genes.use,],2,median)
  cluster_order[k] = which.max(meanA)
}
cluster_order = order(cluster_order)

genes.order = c()
for (c in cluster_order){
  genes.order = c(genes.order, names(a$cluster)[a$cluster==c])
}


p=DataHeatmap(data=t(A[genes.order,]),quantile.thresh = c(0.1,0.9), return.plot = TRUE, cex.row = 0.1, col.low = "blue", col.mid = "white",col.high = "red") 

k=0
for (c in rev(cluster_order)){
  print(c)
  k = k + sum(a$cluster == c)
  p = p + geom_hline(yintercept = k, color="black")
}
p

dir.create("Figs")
pdf("Figs/TemporalHeatmap_RGCdev2.pdf",w=6,h=8, useDingbats = FALSE)
print(p)
dev.off()


cluster_ids = a$cluster
cluster_ids = mapvalues(a$cluster, from = cluster_order, to = c(1:6))
GOres_devRGC = list()
l=1
for (c1 in c(1:6)){
  print(c1)
  res = list()
  
  for (term in c("BP", "MF","CC")){
    print(term)
    res[[term]] = topGOterms(fg.genes = names(cluster_ids)[cluster_ids==c1], bg.genes = rownames(rgc_dev@data), topnodes.print = 800, ontology.use = term)
  }
  
  GOres_devRGC[[c1]] = res
}

dir.create("Txt")
mod_list = list()
for (c1 in c(1:6)){
  mod_list[[c1]] = names(cluster_ids)[cluster_ids == c1]
}

mod_list = lapply(mod_list, `length<-`, max(lengths(mod_list)))
mod_list = lapply(mod_list, function(x) x = c(x[!is.na(x)],rep("", length(x[is.na(x)]))))
names(mod_list) = paste0("Module",c(1:6))
df_module = do.call(cbind.data.frame, mod_list)
rownames(df_module) = NULL
write.table(df_module, file = "Txt/MaturationModulesRGC.txt", sep = "\t", quote=FALSE)
save(list=c("cluster_ids","GOres_devRGC"), file="GOterms_devRGC.Rdata")

# GO heatmaps
for (term in c("BP","MF","CC")){
  
  for (c1 in c(1:6)){
    if (c1 == 1){
      df = GOres_devRGC[[c1]][[term]]$res.table
      df$module = c1
    } else {
      df_temp = GOres_devRGC[[c1]][[term]]$res.table; df_temp$module = c1
      df = rbind(df, df_temp)
    }
  }
  
  # Find terms to plot
  GOterms_to_plot = c()
  for (c1 in c(1:6)){
    tmp_terms = topNonRedundantTerms(GOres_devRGC[[c1]][[term]]$GOdata, 
                                     GOres_devRGC[[c1]][[term]]$res.table,avoid_terms = GOterms_to_plot)
    GOterms_to_plot = unique(c(GOterms_to_plot, tmp_terms))
  }
  
  
  GO_df = -log10(matrix(runif(length(GOterms_to_plot)*6,min=0.1, max=1), nrow=length(GOterms_to_plot)))
  colnames(GO_df) = paste0("Mod",c(1:6))
  rownames(GO_df) = GOterms_to_plot
  l=1
  
  
  row_names_new = c()
  for (c1 in c(1:6)){
    df_temp = subset(df, (GO.ID %in% GOterms_to_plot) & module == c1)
    rownames(df_temp) = df_temp$Term
    pval_log = -log10(as.numeric(df_temp$pval))
    pval_log[is.na(pval_log)] = 30
    GO_df[df_temp$GO.ID, c1] = pval_log
  }
  
  GO_df[GO_df > 5] = 5
  
  GO_terms = rownames(GO_df)
  
  rownames_df = subset(df, GO.ID %in% rownames(GO_df))[,c(1,2)]
  rownames_df = rownames_df[!duplicated(rownames_df),]
  rownames(rownames_df) = rownames_df$GO.ID
  rownames(GO_df) = rownames_df[rownames(GO_df), "Term"]
  pdf(paste0("Figs/", term, "_GOheatmap.pdf"),w=8,h=6,useDingbats = FALSE)
  DataHeatmap(t(GO_df))
  dev.off()
  
  
  # Print table
  gene_terms_all = genesInTerm(GOres_devRGC[[c1]][[term]]$GOdata)
  names(GO_terms) = rownames(GO_df)
  l=1
  for (goterm in GO_terms){
    df_print_tmp = data.frame(GO.ID = goterm, Term = names(GO_terms)[GO_terms==goterm])
    rownames(df_print_tmp) = NULL
    for (c1 in c(1:6)){
      genes_in_mod = intersect(mod_list[[c1]], gene_terms_all[[goterm]])
      df_print_tmp[,paste0("Module",c1)] = paste(genes_in_mod, collapse=",")
    }
    
    if (l==1){
      df_GOprint = df_print_tmp
      l=2
    } else {
      df_GOprint = rbind(df_GOprint, df_print_tmp)
    }
    
  }
  
  write.table(df_GOprint, file=paste0("Txt/", term, "_MaturationGOterms.txt"), sep="\t", quote=FALSE)
  
  
}

topNonRedundantTerms = function(GO_object, 
                                res_table, 
                                topN = 6, 
                                avoid_terms = NULL, 
                                maxN = 200,
                                minN = 5,
                                pval_thresh = 1e-2){
 
  
  gene_terms_all = genesInTerm(GO_object)
  GO_terms_all = c()
  countn = 0
  
  for (idx in c(1:nrow(res_table))){
    GO_tmp = res_table[idx,"GO.ID"]
    Nannot = res_table[idx,"Annotated"]
    pval = as.numeric(res_table[idx,"pval"])
    if (is.na(pval)) pval = 1e-30
    if (pval > pval_thresh) break
    
    if ((Nannot >= maxN) | (Nannot <= minN) ) next
    if ((GO_tmp %in% avoid_terms)) next
    
    if (length(GO_terms_all) == 0){
      GO_terms_all = c(GO_terms_all, GO_tmp)
    } else {
      
      genes_in_term0 = gene_terms_all[[GO_tmp]]
      
      flag_tmp = 0
      for (gterm in GO_terms_all){
        genes_in_term1 = gene_terms_all[[gterm]]
        
        if ((length(intersect(genes_in_term0, genes_in_term1)) / min(length(genes_in_term0),length(genes_in_term1))) > 0.3 ){
          flag_tmp=1
        }
          
      }
      if (flag_tmp == 0) GO_terms_all = c(GO_terms_all, GO_tmp)
      
    }
    
    if (length(GO_terms_all) >= topN) break
    
  }
  
  return(GO_terms_all)
  
}

