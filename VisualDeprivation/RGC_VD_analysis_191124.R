assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")

load("Objects/rgc_perturb_finalfilt_191122.Rdata")
#Read atlas
rgc_atlas = readRDS("../Developing_RGCs/Objects/rgc_atlas_190911_withP56.rds")


# ARI and conditional entropy
order1 = sample(1:length(rgc_perturb@cell.names))
infotheo::condentropy(rgc_perturb@data.info$m_new,rgc_perturb@data.info$iGraphBoostC) / infotheo::entropy(rgc_perturb@data.info$m_new)
mclust::adjustedRandIndex(rgc_perturb@data.info$m_new, rgc_perturb@data.info$iGraphBoostC)
mclust::adjustedRandIndex(rgc_perturb@data.info$m_new[order1], rgc_perturb@data.info$iGraphBoostC)


A = table(rgc_perturb@data.info$iGraphBoostC, rgc_perturb@data.info$m_new)
A = t(scale(t(A), scale=rowSums(A), center=FALSE))
test <- heatmap.2(A, hclustfun = function(x) hclust(x,method="single"))
B=A[test$rowInd, test$colInd]
row.max = apply(B,1,which.max)
row.ord = names(sort(row.max))
pdf("Figs/VD_ConfusionMatrix.pdf", w=8,h=7, useDingbats = FALSE)
plotConfusionMatrix(B[row.ord,], col.low = "white", col.high = "darkblue", ylab.use = "Atlas types", xlab.use = "VD Clusters")
dev.off()




# Relative frequencies
df = data.frame(atlas = table(rgc_atlas@data.info$new_types))
df$RD1 = as.numeric(table(rgc_perturb@data.info$iGraphBoostC[rgc_perturb@data.info$batch == "RD1"])[-46])
df$DR = as.numeric(table(rgc_perturb@data.info$iGraphBoostC[rgc_perturb@data.info$batch == "DR"])[-46])
df$BCless = as.numeric(table(rgc_perturb@data.info$iGraphBoostC[rgc_perturb@data.info$batch == "BCless"])[-46])
rownames(df) = df$atlas.Var1
df = df[,-1]
colnames(df)[1] = "atlas"
df$atlas_b1 = as.numeric(table(rgc_atlas@data.info$new_types[rgc_atlas@data.info$batch == "Batch2"]))
df$atlas_b2 = as.numeric(table(rgc_atlas@data.info$new_types[rgc_atlas@data.info$batch == "Batch3"]))

df$atlas = as.numeric(table(rgc_atlas@data.info$new_types))

df = scale(df, center=FALSE, scale=colSums(df))

df = as.data.frame(df)

library(philentropy)

pdf("Figs/JS_atlas_rd1_freq.pdf",w=6,h=5, useDingbats = FALSE)
ggplot(df, aes(x=atlas, y=RD1)) + geom_ribbon(aes(ymin = 0.33*atlas, ymax = 3*atlas), fill="red", alpha=0.3) + geom_ribbon(aes(ymin = 0.5*atlas, ymax = 2*atlas), fill="yellow", alpha=0.3) + geom_point() + geom_abline(slope = 1, color="gray40", linetype="dashed") + scale_x_log10() + scale_y_log10() + ggtitle(paste0("JSD = ", round(JSD(rbind(df[,"atlas"], df[,"RD1"])),2)))  + theme_classic()
dev.off()

pdf("Figs/JS_atlas_dr_freq.pdf",w=6,h=5, useDingbats = FALSE)
ggplot(df, aes(x=atlas, y=DR)) + geom_ribbon(aes(ymin = 0.33*atlas, ymax = 3*atlas), fill="red", alpha=0.3) + geom_ribbon(aes(ymin = 0.5*atlas, ymax = 2*atlas), fill="yellow", alpha=0.3) + geom_point() + geom_abline(slope = 1, color="gray40", linetype="dashed") + scale_x_log10() + scale_y_log10() + ggtitle(paste0("JSD = ", round(JSD(rbind(df[,"atlas"], df[,"DR"])),2)))  + theme_classic()
dev.off()

pdf("Figs/JS_atlas_bcless_freq.pdf",w=6,h=5, useDingbats = FALSE)
ggplot(df, aes(x=atlas, y=BCless)) + geom_ribbon(aes(ymin = 0.33*atlas, ymax = 3*atlas), fill="red", alpha=0.3) + geom_ribbon(aes(ymin = 0.5*atlas, ymax = 2*atlas), fill="yellow", alpha=0.3) + geom_point() + geom_abline(slope = 1, color="gray40", linetype="dashed") + scale_x_log10() + scale_y_log10() + ggtitle(paste0("JSD = ", round(JSD(rbind(df[,"atlas"], df[,"BCless"])),2)))  + theme_classic()
dev.off()

pdf("Figs/JS_atlas_b1_b2.pdf",w=6,h=5, useDingbats = FALSE)
ggplot(df, aes(x=atlas_b1, y=atlas_b2)) + geom_ribbon(aes(ymin = 0.33*atlas_b1, ymax = 3*atlas_b1), fill="red", alpha=0.3) + geom_ribbon(aes(ymin = 0.5*atlas_b1, ymax = 2*atlas_b1), fill="yellow", alpha=0.3) + geom_point() + geom_abline(slope = 1, color="gray40", linetype="dashed") + scale_x_log10() + scale_y_log10() + ggtitle(paste0("JSD = ", round(JSD(rbind(df[,"atlas_b1"], df[,"atlas_b2"])),2)))  + theme_classic()
dev.off()


###################################################################################
#### Global Gene Expression Changes with respect to Adult 
##################################################################################
# Cleaned perturb object

# Make a small object
rgc_perturb = set.all.ident(rgc_perturb, id="batch")
cells.use = c(sample(which.cells(rgc_perturb,"RD1"),8000), sample(which.cells(rgc_perturb,"DR"),8000),
              sample(which.cells(rgc_perturb,"BCless"),8000))
cells.use_adult = sample(rgc_atlas@cell.names, 8000)

genes.use = intersect(rownames(rgc_perturb@count.data), rownames(rgc_atlas@count.data))
Count.mat = cbind(rgc_perturb@count.data[genes.use, cells.use], rgc_atlas@count.data[genes.use,cells.use_adult])

rgc_adult_pert=scR(count.data=Count.mat,ident.fxn=getStat1)
rgc_adult_pert=setup(rgc_adult_pert,project="RGC_adultpert",min.cells = 1,min.genes = 500,is.expr=0, threshold.quantile = 0.999)

rgc_adult_pert@data.info$Condition = "Adult"
rgc_adult_pert@data.info[grep("darkrear", rgc_adult_pert@cell.names, value=TRUE), "Condition"] = "DR"
rgc_adult_pert@data.info[grep("RD1", rgc_adult_pert@cell.names, value=TRUE), "Condition"] = "RD1"
rgc_adult_pert@data.info[grep("RGCBCless", rgc_adult_pert@cell.names, value=TRUE), "Condition"] = "BCless"
rgc_adult_pert@data.info$Condition = factor(rgc_adult_pert@data.info$Condition, levels = c("Adult","DR","RD1",
                                                                                           "BCless"))

rgc_adult_pert = set.all.ident(rgc_adult_pert, id="Condition")

DE_vs_adult_global = list()

for (l in c("RD1","DR","BCless")){
  print(l)
  x = find.markers(rgc_adult_pert, l, "Adult", thresh.use = 0.2, test.use="MAST")
  DE_vs_adult_global[[l]] = x
}

save(list=c("DE_vs_adult_global"), file="DE_Perturbation_vs_Adult_191124.Rdata")


rgc_perturb = set.all.ident(rgc_perturb, id="iGraphBoostC")
rgc_atlas@data.info$clustC = rgc_atlas@data.info$new_types
levels(rgc_atlas@data.info$clustC) = paste0("C",c(1:45))

# type specific markers
DE_vs_adult_type = list()
DE_vs_adult_type[["RD1"]]  = list()
DE_vs_adult_type[["DR"]]  = list()
DE_vs_adult_type[["BCless"]]  = list()

pdf("Figs/VlnPlot_nGene_Rbpms.pdf",w=6,h=4)
VlnPlot(rgc_perturb,"nGene")
VlnPlot(rgc_perturb,"nTranscripts")
VlnPlot(rgc_perturb,"Rbpms")
VlnPlot(rgc_perturb,"Pou4f1")
VlnPlot(rgc_perturb,"Slc17a6")
VlnPlot(rgc_perturb,"Thy1")
dev.off()

rgc_atlas = set.all.ident(rgc_atlas, id="clustC")
for (l in levels(rgc_atlas@ident)){
  print(l)
  adult_type = subsetData(rgc_atlas, cells.use = which.cells(rgc_atlas, l))
  pert_type_all = subsetData(rgc_perturb, cells.use = which.cells(rgc_perturb, l))
  pert_type_all = set.all.ident(pert_type_all, id="batch")
  
  for (pert in c("RD1","DR","BCless")){
    pert_type = subsetData(pert_type_all, cells.use = which.cells(pert_type_all, pert))
    DE_vs_adult_type[[pert]][[l]]$atlas_ncell = length(adult_type@cell.names)
    DE_vs_adult_type[[pert]][[l]]$ncell = length(pert_type@cell.names)
    pert_type = set.all.ident(pert_type, id="iGraphBoostC")
    if (length(pert_type@cell.names) < 50){
      DE_vs_adult_type[[pert]][[l]]$DE = NA
      next
    }
    combine_clust = combineObjects(pert_type, adult_type, tag1 = "Pert", tag2 = "Atlas", use.gene.union = TRUE)
    a=find.markers(combine_clust, paste0("Pert_", l),paste0("Atlas_",l), test.use="MAST", thresh.use = 0.005, min.perc = 0.1, min.diff = 0.03)
    DE_vs_adult_type[[pert]][[l]]$DE = a
  }
  save(list=c("DE_vs_adult_global", "DE_vs_adult_type"), file="DE_Perturbation_vs_Adult_191124.Rdata")
  
}

save(list=c("DE_vs_adult_global", "DE_vs_adult_type"), file="DE_Perturbation_vs_Adult_191124.Rdata")
load("DE_Perturbation_vs_Adult_191124.Rdata")

#########################################################
### Analyze Global Differences #########
########################################################
load("DE_Perturbation_vs_Adult_191124.Rdata")

library(UpSetR)
up_genes = list(); down_genes = list()
for (g in c("RD1", "DR", "BCless")){
  up_genes[[g]] = rownames(subset(DE_vs_adult_global[[g]], (diff >= log(1.5)) & (pct_clust1 > 0.7)))
  down_genes[[g]] = rownames(subset(DE_vs_adult_global[[g]], (diff <= -log(1.5)) & (pct_clust2 > 0.7)))
  
  if (g == "RD1"){
    df = subset(DE_vs_adult_global[[g]], abs(diff) >= log(1.25))
    df$condition = g
    df$gene = rownames(df)
    colnames(df)[7] = "nTrans_condition"
  } else {
    df_temp = subset(DE_vs_adult_global[[g]], abs(diff) >= log(1.25))
    df_temp$condition = g
    df_temp$gene = rownames(df_temp)
    colnames(df_temp)[7] = "nTrans_condition"
    df = rbind(df, df_temp)
    
  }
}
rownames(df) = NULL
write.table(df, file = "Objects/DE_global_RD1_DR_BCless_191124.txt", sep="\t", quote=FALSE)

all_up_genes = unique(unlist(up_genes))
all_down_genes = unique(unlist(down_genes))

df_up = as.data.frame(matrix(0, nrow=length(all_up_genes), ncol=3, dimnames = list(all_up_genes, c("RD1", "DR", "BCless"))))
df_down = as.data.frame(matrix(0, nrow=length(all_down_genes), ncol=3, dimnames = list(all_down_genes, c("RD1", "DR", "BCless"))))

for (g in c("RD1", "DR", "BCless")){
  df_up[up_genes[[g]], g] = 1
  df_down[down_genes[[g]], g] = 1
}

pdf("Figs/Global_genes_DR_RD1_BCless_barplot_upset.pdf",w=9,h=5, useDingbats = FALSE)
upset(df_down, point.size = 3)
upset(df_up, point.size = 3)
dev.off()

df_bar = data.frame(n=c(84,96,113,28,10,0,0), sign = "Down", set = c("DR","RD1","BCless",
                                                                       "DR_RD1","RD1_BCless", "DR_BCless",
                                                                       "DR_RD1_BCless"))
df_bar2 = data.frame(n=c(12,10,86,7,11,11,9), sign = "Up", set = c("DR","RD1","BCless",
                                                                      "DR_RD1","RD1_BCless", "DR_BCless",
                                                                      "DR_RD1_BCless"))

df = rbind(df_bar, df_bar2)
df$set = factor(df$set, levels = c("DR","RD1","BCless",
                                   "DR_RD1","RD1_BCless", "DR_BCless",
                                   "DR_RD1_BCless"))

pdf("Figs/Global_genes_DR_RD1_BCless_barplot3.pdf",w=9,h=5, useDingbats = FALSE)
ggplot(data=df, aes(x=set, y=n)) + geom_bar(position="dodge",stat="identity", aes(fill=factor(sign)), color="black") + scale_fill_manual(values=c("red","blue")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("# genes (fold > 1.5)") + theme_classic()
dev.off()



# Look at global differences  (Down)
downGO = list()
downGO[["DR"]] = topGOterms(fg.genes = down_genes[["DR"]], bg.genes = c(all_up_genes, all_down_genes), ontology.use = "BP")
downGO[["RD1"]] = topGOterms(fg.genes = down_genes[["RD1"]], bg.genes = c(all_up_genes, all_down_genes), ontology.use = "BP")

downGO[["BCless"]] = topGOterms(fg.genes = down_genes[["BCless"]], bg.genes = c(all_up_genes, all_down_genes), ontology.use = "BP")

df = downGO[["DR"]]$res.table; df$condition = "DR"
df_temp = downGO[["RD1"]]$res.table; df_temp$condition = "RD1"; df = rbind(df, df_temp)
df_temp = downGO[["BCless"]]$res.table; df_temp$condition = "BCless";df = rbind(df, df_temp)

terms = downGO[["DR"]]$res.table$GO.ID[c(1,2,3,5,9)]
terms = c(terms,downGO[["RD1"]]$res.table$GO.ID[c(1,2,4,5)] )
terms = c(terms,downGO[["BCless"]]$res.table$GO.ID[c(2,3,7)] )

downGO_VD = -log10(matrix(runif(12*3,min=0.1, max=1), nrow=12))
colnames(downGO_VD) = c("DR","RD1","BCless")
rownames(downGO_VD) = terms
l=1

row_names_new = c()
for (g in c("DR","RD1","BCless")){
  df_temp = subset(df, (GO.ID %in% terms) & condition == g)
  rownames(df_temp) = df_temp$Term
  row_names_new = c(row_names_new, rownames(df_temp))
  downGO_VD[df_temp$GO.ID, g] = -log10(as.numeric(df_temp$pval))
}
pdf("Figs/GlobalGO_VD_down.pdf",w=7,h=6,useDingbats = FALSE)
DataHeatmap(t(DownGO_VD))
dev.off()

# Up genes

upGO = list()
upGO[["DR"]] = topGOterms(fg.genes = up_genes[["DR"]], bg.genes = c(all_up_genes, all_down_genes), ontology.use = "BP")
upGO[["RD1"]] = topGOterms(fg.genes = up_genes[["RD1"]], bg.genes = c(all_up_genes, all_down_genes), ontology.use = "BP")

upGO[["BCless"]] = topGOterms(fg.genes = up_genes[["BCless"]], bg.genes = c(all_up_genes, all_down_genes), ontology.use = "BP")

df = upGO[["DR"]]$res.table; df$condition = "DR"
df_temp = upGO[["RD1"]]$res.table; df_temp$condition = "RD1"; df = rbind(df, df_temp)
df_temp = upGO[["BCless"]]$res.table; df_temp$condition = "BCless";df = rbind(df, df_temp)

terms = upGO[["DR"]]$res.table$Term[c(1,2,3,6,7,9)]
terms = c(terms,upGO[["RD1"]]$res.table$Term[c(1,3,4,5,6,8,9)] )
terms = c(terms,upGO[["BCless"]]$res.table$Term[c(1,2,3,5,14)] )
terms = unique(terms)



upGO_VD = -log10(matrix(runif(15*3,min=0.1, max=1), nrow=15))
colnames(upGO_VD) = c("DR","RD1","BCless")
rownames(upGO_VD) = terms
l=1

for (g in c("DR","RD1","BCless")){
  df_temp = subset(df, (Term %in% terms) & condition == g)
  rownames(df_temp) = df_temp$Term
  upGO_VD[df_temp$Term, g] = -log10(as.numeric(df_temp$pval))
}

pdf("Figs/GlobalGO_VD_up.pdf",w=6.5,h=6,useDingbats = FALSE)
DataHeatmap(t(upGO_VD), quantile.thresh = c(0.1, 0.83))
dev.off()



#########################################################
### Analyze Type Specific Differences #########
########################################################
load("DE_Perturbation_vs_Adult_191124.Rdata")

library(UpSetR)
# DE genes vs. type size
rgc_types = names(DE_vs_adult_type$RD1)
cond = c("RD1", "DR", "BCless")
df = data.frame(type = as.character(), condition = as.character(), num_genes = as.numeric(), num_cells = as.numeric())

l=1
for (c1 in cond){
  de_global = subset(DE_vs_adult_global[[c1]], abs(diff) >= log(1.5) & (pct_clust1 > 0.7 | pct_clust2 > 0.7) )
  for (rgc_type in rgc_types){
    de = DE_vs_adult_type[[c1]][[rgc_type]]
    if (is.na(de$DE)) next
    de$DE = de$DE[setdiff(rownames(de$DE), intersect(rownames(de$DE), rownames(de_global))),]
    df = rbind(df,data.frame(type = rgc_type, condition = c1, num_genes = nrow(subset(de$DE, abs(diff) >= log(1.5) )),num_cells = min(de$atlas_ncell, de$ncell)))
    
    if (l==1){
      df_temp = subset(de$DE, abs(diff) >= log(1.2))
      colnames(df_temp)[7] = "nTrans_condition"
      colnames(df_temp)[8] = "nTrans_Atlas"
      df_temp$type = rgc_type
      df_temp$condition = c1
      df_temp$gene = rownames(df_temp)
      df_print = df_temp
      l=l+1
    } else {
      df_temp = subset(de$DE, abs(diff) >= log(1.2))
      colnames(df_temp)[7] = "nTrans_condition"
      colnames(df_temp)[8] = "nTrans_Atlas"
      df_temp$type = rgc_type
      df_temp$condition = c1
      df_temp$gene = rownames(df_temp)
      df_print = rbind(df_print, df_temp)
    }
  }
}

rownames(df_print) = NULL
write.table(df_print, file = "Objects/DE_type_RD1_DR_BCless_191124.txt", sep="\t", quote=FALSE)

pdf("Figs/Num_DE_genes_vs_cells_for_type2.pdf",w=6,h=4.5, useDingbats = FALSE)
ggplot(subset(df, num_cells < 6000), aes(x=num_cells, y=num_genes)) + geom_point(aes(color=condition)) +  geom_smooth() + xlab("Cells per type") + ylab("# DE genes (log-fold > 1.5)") + scale_x_log10() + scale_y_log10() + theme_classic()
dev.off()

# Overlap of type specific genes across perturbations
library(UpSetR)
up_genes = list(); down_genes = list()
for (g in c("RD1", "DR", "BCless")){
  up_genes[[g]] = c()
  down_genes[[g]] = c()
  de_global = subset(DE_vs_adult_global[[g]], abs(diff) >= log(1.5) & (pct_clust1 > 0.7 | pct_clust2 > 0.7))
  
  for (rgc_type in rgc_types){
    de = DE_vs_adult_type[[g]][[rgc_type]]
    if (is.na(de$DE)) next
    de$DE = de$DE[setdiff(rownames(de$DE), intersect(rownames(de$DE), rownames(de_global))),]
    
    up_genes[[g]] = union( up_genes[[g]], rownames(subset(de$DE, diff >= log(1.5) & pct_clust1 > 0.3))) 
    down_genes[[g]] = union( down_genes[[g]], rownames(subset(de$DE, diff <= -log(1.5) & pct_clust2 > 0.3))) 
  }
  
}

all_up_genes = unique(unlist(up_genes))
all_down_genes = unique(unlist(down_genes))

df_up = as.data.frame(matrix(0, nrow=length(all_up_genes), ncol=3, dimnames = list(all_up_genes, c("RD1", "DR", "BCless"))))
df_down = as.data.frame(matrix(0, nrow=length(all_down_genes), ncol=3, dimnames = list(all_down_genes, c("RD1", "DR", "BCless"))))

for (g in c("RD1", "DR", "BCless")){
  df_up[up_genes[[g]], g] = 1
  df_down[down_genes[[g]], g] = 1
}
upset(df_up)
upset(df_down)
df_bar = data.frame(n=c(1314,482,638,267,169,87,14), sign = "Down", set = c("DR","RD1","BCless",
                                                                             "DR_RD1","RD1_BCless", "DR_BCless",
                                                                             "DR_RD1_BCless"))
df_bar2 = data.frame(n=c(154,192,912,32,93,42,18), sign = "Up", set = c("DR","RD1","BCless",
                                                                          "DR_RD1","RD1_BCless", "DR_BCless",
                                                                          "DR_RD1_BCless"))

df = rbind(df_bar, df_bar2)
df$set = factor(df$set, levels = c("DR","RD1","BCless",
                                   "DR_RD1","RD1_BCless", "DR_BCless",
                                   "DR_RD1_BCless"))

pdf("Figs/TypeSp_genes_DR_RD1_BCless_barplot2.pdf",w=9,h=5, useDingbats = FALSE)
ggplot(data=df, aes(x=set, y=n)) + geom_bar(position="dodge",stat="identity", aes(fill=factor(sign)), color="black") + scale_fill_manual(values=c("red","blue")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("# DE genes (fold > 1.5)") + theme_classic()
dev.off()



# Sharing of type specific genes (some code here is repeated from above)

up_genes_type = list(); down_genes_type = list()
up_genes_global = list(); down_genes_global = list()
for (g in c("RD1", "DR", "BCless")){
  up_genes_type[[g]] = c()
  down_genes_type[[g]] = c()
  up_genes_global[[g]] = c()
  down_genes_global[[g]] = c()
  de_global = subset(DE_vs_adult_global[[g]], abs(diff) >= log(1.5) & (pct_clust1 > 0.7 | pct_clust2 > 0.7))
  global_genes = rownames(de_global)
  
  for (rgc_type in rgc_types){
    de = DE_vs_adult_type[[g]][[rgc_type]]
    if (is.na(de$DE)) next
    #de$DE = de$DE[setdiff(rownames(de$DE), intersect(rownames(de$DE), rownames(de_global))),]
    up_genes_type[[g]] = c( up_genes_type[[g]], rownames(subset(de$DE, diff >= log(1.5) & pct_clust1 > 0.3))) 
    down_genes_type[[g]] = c( down_genes_type[[g]], rownames(subset(de$DE, diff <= -log(1.5) & pct_clust2 > 0.3))) 
  }
  
  up_genes_global[[g]] = up_genes_type[[g]][up_genes_type[[g]] %in% global_genes]
  up_genes_type[[g]] = up_genes_type[[g]][!(up_genes_type[[g]] %in% global_genes)]
  down_genes_global[[g]] = down_genes_type[[g]][down_genes_type[[g]] %in% global_genes]
  down_genes_type[[g]] = down_genes_type[[g]][!(down_genes_type[[g]] %in% global_genes)]
  
}

for (g in c("RD1", "DR", "BCless")){
  
  # Type specific Up
  df_temp = data.frame(x=c(1:45), n=0, condition = g, direction = "Up", type="specific")
  gene_freq = table(up_genes_type[[g]])
  for (ind in c(1:45)){
    df_temp[df_temp$x==ind,"n"] = sum(gene_freq==ind)
  }
  
  a = df_temp[,"n"]
  df_temp[,"n"] = cumsum(a) / sum(a)
  
  # Type specific Down
  df_temp2 = data.frame(x=c(1:45), n=0, condition = g, direction = "Down", type="specific")
  gene_freq = table(down_genes_type[[g]])
  for (ind in c(1:45)){
    df_temp2[df_temp2$x==ind,"n"] = sum(gene_freq==ind)
  }
  a = df_temp2[,"n"]
  df_temp2[,"n"] = cumsum(a) / sum(a)
  
  # Global Up
  df_temp3 = data.frame(x=c(1:45), n=0, condition = g, direction = "Up", type="global")
  gene_freq = table(up_genes_global[[g]])
  for (ind in c(1:45)){
    df_temp3[df_temp$x==ind,"n"] = sum(gene_freq==ind)
  }
  
  a = df_temp3[,"n"]
  df_temp3[,"n"] = cumsum(a) / sum(a)
  
  # Global down
  df_temp4 = data.frame(x=c(1:45), n=0, condition = g, direction = "Down", type="global")
  gene_freq = table(down_genes_global[[g]])
  for (ind in c(1:45)){
    df_temp4[df_temp2$x==ind,"n"] = sum(gene_freq==ind)
  }
  a = df_temp4[,"n"]
  df_temp4[,"n"] = cumsum(a) / sum(a)
  
  if (g=="RD1"){ 
    df = rbind(rbind(df_temp, df_temp2), rbind(df_temp3, df_temp4))
  } else {
    df = rbind(df, rbind(rbind(df_temp, df_temp2), rbind(df_temp3, df_temp4)))
  }
  
}

df$fac = paste0(as.character(df$condition),"_", as.character(df$direction),"_", as.character(df$type))
l1 = paste0(c("DR_Up","RD1_Up","BCless_Up","DR_Down","RD1_Down","BCless_Down"),"_specific")
l2 = paste0(c("DR_Up","RD1_Up","BCless_Up","DR_Down","RD1_Down","BCless_Down"),"_global")

df$fac = factor(df$fac, levels = c(l1,l2))
pdf("Figs/TypesPerDEgene_CDF.pdf",w=6,h=5, useDingbats = FALSE)
ggplot(df, aes(x=x, y=n, group=fac)) + geom_point(aes(color=direction, shape=condition)) + geom_line(aes(color=direction, linetype=type)) + scale_color_manual(values = c("red","blue")) + theme_classic() + xlab("# types per gene") + ylab("CDF(x)")
dev.off()



topGOterms(fg.genes = rownames(subset(df_down, RD1==1)), bg.genes = c(unique(rownames(df_up), rownames(df_down))), ontology.use = "BP")

#############
### GO analysis
#############

all_up_genes_type = unique(unlist(up_genes_type))
all_up_genes_global = unique(unlist(up_genes_global))
all_down_genes_type = unique(unlist(down_genes_type))
all_down_genes_global = unique(unlist(down_genes_global))

# Look at global differences  (Down)
bg_genes = c(all_up_genes_type, all_down_genes_type, all_up_genes_global, all_down_genes_global)
downGO = list()
downGO[["DR"]] = topGOterms(fg.genes = c(down_genes_global[["DR"]],down_genes_type[["DR"]]) , bg.genes = bg_genes, ontology.use = "BP", topnodes.print = 500)
downGO[["RD1"]] = topGOterms(fg.genes = c(down_genes_global[["RD1"]],down_genes_type[["RD1"]]), bg.genes = bg_genes, ontology.use = "BP",topnodes.print = 500)

downGO[["BCless"]] = topGOterms(fg.genes = c(down_genes_global[["BCless"]],down_genes_type[["BCless"]]), bg.genes = bg_genes, ontology.use = "BP",topnodes.print = 500)

df = downGO[["DR"]]$res.table; df$condition = "DR"
df_temp = downGO[["RD1"]]$res.table; df_temp$condition = "RD1"; df = rbind(df, df_temp)
df_temp = downGO[["BCless"]]$res.table; df_temp$condition = "BCless";df = rbind(df, df_temp)

terms = downGO[["DR"]]$res.table$Term[c(2,3,4,6,7,16,19)]
terms = c(terms,downGO[["RD1"]]$res.table$Term[c(4,5,6,10)] )
terms = c(terms,downGO[["BCless"]]$res.table$Term[c(2,3,8,9,10)] )
terms = unique(terms)
downGO_VD = -log10(matrix(runif(length(terms)*3,min=0.1, max=1), nrow=length(terms)))
colnames(downGO_VD) = c("DR","RD1","BCless")
rownames(downGO_VD) = terms
l=1

row_names_new = c()
for (g in c("DR","RD1","BCless")){
  df_temp = subset(df, (Term %in% terms) & condition == g)
  rownames(df_temp) = df_temp$Term
  row_names_new = c(row_names_new, rownames(df_temp))
  downGO_VD[df_temp$Term, g] = -log10(as.numeric(df_temp$pval))
}

pdf("Figs/GlobalTypeGO_VD_down.pdf",w=7,h=6,useDingbats = FALSE)
DataHeatmap(t(downGO_VD))
dev.off()

# Up genes

upGO = list()
upGO[["DR"]] = topGOterms(fg.genes = c(up_genes_type[["DR"]], up_genes_global[["DR"]]), bg.genes = bg_genes, ontology.use = "BP", topnodes.print = 500)
upGO[["RD1"]] = topGOterms(fg.genes = c(up_genes_type[["RD1"]], up_genes_global[["RD1"]]), bg.genes = bg_genes, ontology.use = "BP", topnodes.print = 500)

upGO[["BCless"]] = topGOterms(fg.genes = c(up_genes_type[["BCless"]], up_genes_global[["BCless"]]), bg.genes = bg_genes, ontology.use = "BP", topnodes.print = 500)

df = upGO[["DR"]]$res.table; df$condition = "DR"
df_temp = upGO[["RD1"]]$res.table; df_temp$condition = "RD1"; df = rbind(df, df_temp)
df_temp = upGO[["BCless"]]$res.table; df_temp$condition = "BCless";df = rbind(df, df_temp)

terms = upGO[["DR"]]$res.table$Term[c(2,3,4,6,7,9)]
terms = c(terms,upGO[["RD1"]]$res.table$Term[c(2,5,7,9,12)] )
terms = c(terms,upGO[["BCless"]]$res.table$Term[c(2,3,4,5,10,11)] )
terms = unique(terms)



upGO_VD = -log10(matrix(runif(length(terms)*3,min=0.1, max=1), nrow=length(terms)))
colnames(upGO_VD) = c("DR","RD1","BCless")
rownames(upGO_VD) = terms
l=1

for (g in c("DR","RD1","BCless")){
  df_temp = subset(df, (Term %in% terms) & condition == g)
  rownames(df_temp) = df_temp$Term
  upGO_VD[df_temp$Term, g] = -log10(as.numeric(df_temp$pval))
}

pdf("Figs/GlobalTypeGO_VD_up.pdf",w=6.5,h=6,useDingbats = FALSE)
DataHeatmap(t(upGO_VD))
dev.off()


