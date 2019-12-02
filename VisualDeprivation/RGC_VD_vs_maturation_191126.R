assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")


#### Atlas
rgc_atlas = readRDS("../Developing_RGCs/Objects/rgc_atlas_190911_withP56.rds")

#### P5
rgcP5 = readRDS("../Developing_RGCs/Objects/rgcP5_190911_withP56.rds")

# Load P5 rheaume object
load("../Developing_RGCs/Objects/p5_Rheaumme_object.Rdata")

# RGC perturbation
load("Objects/rgc_perturb_finalfilt_191122.Rdata")

###### Consistency between P5 and Rheaumme
# First make a combined object of the two P5s and the adult
Count.mat_p5rhe = p5rhe@count.data; colnames(Count.mat_p5rhe) = paste0("P5Rhe_", colnames(Count.mat_p5rhe))
Count.mat_p5 = rgcP5@count.data; colnames(Count.mat_p5) = paste0("P5_", colnames(Count.mat_p5))
Count.mat_adult = rgc_atlas@count.data; colnames(Count.mat_adult) = paste0("Adult_", colnames(Count.mat_adult))
genes.use = intersect(rownames(Count.mat_p5rhe), rownames(Count.mat_p5))
genes.use = intersect(genes.use, rownames(Count.mat_adult))

Count.mat = cbind(Count.mat_adult[genes.use, ], Count.mat_p5[genes.use,])
Count.mat = cbind(Count.mat, Count.mat_p5rhe[genes.use,])

p5adult=scR(count.data=Count.mat,ident.fxn=getStat1)
p5adult=setup(p5adult,project="P5_adult",min.cells = 1,min.genes = 500,is.expr=0, threshold.quantile = 0.999)
a=find.markers(p5adult,"Adult", "P5",test.use="MAST", max.cells.per.ident = 3000)
b=find.markers(p5adult,"Adult", "P5Rhe",test.use="MAST", max.cells.per.ident = 3000)

genes.use = union(rownames(a), rownames(b))
ExpMat = Perc.pos.by.ident(p5adult, rownames(p5adult@data), do.plot = FALSE)
ExpStrength = ExpMat$PercMat * ExpMat$ExpMat
f1 = (ExpStrength[,1] + 0.01) / (ExpStrength[,2] + 0.01)
f2 = (ExpStrength[,1] + 0.01) / (ExpStrength[,3] + 0.01)

dir.create("Figs")
df = data.frame(f1 = log2(f1), f2 = log2(f2))
pdf("Figs/P5_vs_Adult_consistency_with_Rheaume.pdf",w=7,h=6, useDingbats = FALSE)
ggplot(df, aes(f1,f2)) + geom_point(size=0.8, color="gray40") + geom_density_2d(color="orange") + 
  xlab("log2(Expression at Adult / Expression at P5 (This study))") + 
  ylab("log2(Expression at Adult / Expression at P5 (Rheaume et al.))") + 
  ggtitle("R = 0.65") + theme_classic()
dev.off()

df = data.frame(f1 = log2(ExpStrength[,2] + 0.01), f2 = log2(ExpStrength[,3] + 0.01))
pdf("Figs/P5_expression_consistency_with_Rheaume.pdf",w=7,h=6, useDingbats = FALSE)
ggplot(df, aes(f1,f2)) + geom_point(size=0.8, color="gray40") + geom_density_2d(color="orange") + 
  xlab("log2(P5 Expression (This Study))") + 
  ylab("log2(P5 Expression (Rheaume et al.))") + 
  ggtitle("R = 0.91") + theme_classic()
dev.off()

rm(p5adult); rm(p5rhe)

#######################################
# Global changes between datasets (RGC perturbations)
#######################################
Count.mat_adult = rgc_atlas@count.data; colnames(Count.mat_adult) = paste0("Adult_", colnames(Count.mat_adult))
Count.mat_p5 = rgcP5@count.data; colnames(Count.mat_p5) = paste0("P5_", colnames(Count.mat_p5))
Count.mat_pert = rgc_perturb@count.data
colnames(Count.mat_pert) = paste0(as.character(rgc_perturb@data.info$batch),"_", colnames(Count.mat_pert))


genes.use = intersect(rownames(Count.mat_adult), rownames(Count.mat_p5))
genes.use = intersect(genes.use, rownames(Count.mat_pert))

Count.mat = cbind(Count.mat_adult[genes.use, ], Count.mat_p5[genes.use,])
Count.mat = cbind(Count.mat, Count.mat_pert[genes.use,])

rm(Count.mat_adult); rm(Count.mat_pert);  rm(Count.mat_p5);
Count.mat_sample = Count.mat[,sample(colnames(Count.mat),60000)]

rgc_devpert=scR(count.data=Count.mat_sample,ident.fxn=getStat1)
rgc_devpert=setup(rgc_devpert,project="RGC_developmentalpert",min.cells = 1,min.genes = 500,is.expr=0, threshold.quantile = 0.999)

ExpMat = Perc.pos.by.ident(rgc_devpert, rownames(rgc_devpert@data), do.plot = FALSE)
ExpStrength = ExpMat$PercMat * ExpMat$ExpMat

col.order = c("P5", "Adult", "DR", "RD1","BCless")
C=cor(ExpStrength[,col.order])
p=DataHeatmap(data=C, return.plot = TRUE, cex.row = 0.1, col.low = "white", col.mid = "white",col.high = "darkgreen") 

pdf("Figs/CorrelationMatrix_allRGCperturbations.pdf",w=7,h=5,useDingbats = FALSE)
p
dev.off()


# DE genes
bulkDE_VDP5 = list()

for (c1 in c("P5","RD1","DR","BCless")){
  bulkDE_VDP5[[c1]] = list()
  
  bulkDE_VDP5[[c1]]$all = find.markers(rgc_devpert, "Adult",c1,thresh.use = 0.005, test.use="MAST", max.cells.per.ident = 3000)
  bulkDE_VDP5[[c1]]$sig = subset(bulkDE_VDP5[[c1]]$all, (diff > 0 & pct_clust1 > 0.3) | (diff < 0 & pct_clust2 > 0.3))
  
}

save(list=c("bulkDE_VDP5"), file="VD_P5_DE_results.Rdata")

# Scatter plots
for (c1 in c("RD1","DR","BCless")){

  genes.use = union(rownames(bulkDE_VDP5[[c1]]$sig),rownames(bulkDE_VDP5[["P5"]]$sig))
  ExpMat = Perc.pos.by.ident(rgc_devpert, genes.use, do.plot = FALSE)
  ExpStrength = ExpMat$PercMat * ExpMat$ExpMat
  ExpStrength = ExpStrength + 0.01
  
  X1 = sapply(genes.use, function(x){
    if (ExpMat$PercMat[x,"Adult"] < 0.3 & ExpMat$PercMat[x,"P5"] < 0.3){
      return(NA)
    } else {
      return(ExpStrength[x,"Adult"]/ExpStrength[x,"P5"])
    }
  })
  
  X2 = sapply(genes.use, function(x){
    if (ExpMat$PercMat[x,"Adult"] < 0.3 & ExpMat$PercMat[x,c1] < 0.3){
      return(NA)
    } else {
      return(ExpStrength[x,"Adult"]/ExpStrength[x,c1])
    }
  })
  df = data.frame(P5 = log(X1))
  df[,c1] = log(X2)
  df=df[!is.na(rowSums(df)),]
  # Octets
  x1 = sum(df[,1] <= -0.693 & abs(df[,2]) < 0.693)
  x2 = sum(df[,1] <= -0.693 & df[,2] >= 0.693)
  x3 = sum(abs(df[,1]) < 0.693 & df[,2] >= 0.693)
  x4 = sum(df[,1] >= 0.693 & df[,2] >= 0.693)
  x5 = sum(df[,1] >= 0.693 & abs(df[,2]) < 0.693)
  x6 = sum(df[,1] >= 0.693 & df[,2] <= -0.693)
  x7 = sum(abs(df[,1]) < 0.693 & df[,2] <= -0.693)
  x8 = sum(df[,1] < -0.693 & df[,2] < -0.693)
  
  
  pdf(paste0("Figs/",c1,"_vs_P5_DEscatter.pdf"),w=6,h=5, useDingbats = FALSE)
  ggplot(df, aes_string(x="P5",y=c1)) + geom_point(size=0.4, color="gray20") + xlim(c(-5,5)) + ylim(c(-4,4)) + geom_hline(yintercept = 0.693, color="red", linetype="dashed") + 
    geom_hline(yintercept = -0.693, color="red", linetype="dashed") + geom_vline(xintercept = 0.693, color="red", linetype="dashed") + geom_vline(xintercept = -0.693, color="red", linetype="dashed") + 
    xlab("Log-Fold Change (Adult vs P5)") + ylab("Log-Fold Change (Adult vs. RD1)") + 
    annotate("text", x=-4, y=0, label = paste0("n = ", x1), size=5, color="blue") + 
    annotate("text", x= -2, y=2, label = paste0("n = ", x2), size=5, color="blue")  + 
    annotate("text", x= 0, y=2, label = paste0("n = ", x3), size=5, color="blue") +
    annotate("text", x= 2, y=2, label = paste0("n = ", x4), size=5, color="blue")  +
    annotate("text", x= 4, y=0, label = paste0("n = ", x5), size=5, color="blue")  +
    annotate("text", x= 2, y=-2, label = paste0("n = ", x6), size=5, color="blue")  +
    annotate("text", x= 0, y=-2, label = paste0("n = ", x7), size=5, color="blue")  +
    annotate("text", x= -2, y=-2, label = paste0("n = ", x8), size=5, color="blue") + theme_classic()
  dev.off()
}




