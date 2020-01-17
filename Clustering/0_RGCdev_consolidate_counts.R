assign("dirpath", "/home/kshekhar/CompRepos/single_cellR/", envir = .GlobalEnv)
assign("fast_tsne_script", "fast_tsne.R", envir = .GlobalEnv)

source(paste0(dirpath,"ObjectConversion.R"))
source(paste0(dirpath,"sc.R"))
source(paste0(dirpath, "scFunctions.R"))
source(paste0(dirpath,"scR_DE.R"))
source(paste0(dirpath,"scR_plotting.R"))
UMAP <- reticulate::import("umap")


################################################
##### Read Data ###############################
###############################################

# Load Dec 2017 data
load("../../RawData/mouseP0E16_Dec2017_10X.Rdata")
colnames(Count.mat_mouseP0E16) = gsub("RGC_","RGCS1_",colnames(Count.mat_mouseP0E16))

# Load Feb 2018 data
load("../../RawData/MouseEarlyRGC_Feb09_10X.Rdata")
colnames(Count.mat_EarlyRGC) = gsub("RGC_","RGCS1_",colnames(Count.mat_EarlyRGC))

Count.mat = cbind(Count.mat_mouseP0E16, Count.mat_EarlyRGC)
rm(Count.mat_EarlyRGC); rm(Count.mat_mouseP0E16)

# Rename P0S2 as E14 s2
cells.p0_1 = grep("P0Cd90RGCS2", colnames(Count.mat))
cells.p0_2 = grep("P0L1camRGCS2", colnames(Count.mat))

cells.e14_1 = grep("E14Cd90RGCS1", colnames(Count.mat))
cells.e14_2 = grep("E14L1camRGCS1", colnames(Count.mat))

colnames(Count.mat)[cells.p0_1] = gsub("P0Cd90RGCS2","E14Cd90RGCS1",colnames(Count.mat)[cells.p0_1])
colnames(Count.mat)[cells.p0_2] = gsub("P0L1camRGCS2","E14L1camRGCS1",colnames(Count.mat)[cells.p0_2])
colnames(Count.mat)[cells.e14_1] = gsub("E14Cd90RGCS1","P0Cd90RGCS2",colnames(Count.mat)[cells.e14_1])
colnames(Count.mat)[cells.e14_2] = gsub("E14L1camRGCS1","P0L1camRGCS2",colnames(Count.mat)[cells.e14_2])

# Add new E13E14 data
NewData = readRDS("../../RawData/RGC_E13E14_July2019.rds")
E14_new = NewData$Counts_e14s2; E13_new = NewData$Counts_e13s1
colnames(E14_new) = gsub("E14CD9021","E14Cd90RGCS2", colnames(E14_new))
colnames(E14_new) = gsub("L1camS2","L1camRGCS2", colnames(E14_new))
colnames(E13_new) = gsub("CD90S1","Cd90RGCS1", colnames(E13_new))
colnames(E13_new) = gsub("L1camS1","L1camRGCS1", colnames(E13_new))

# Collapse duplicate rows
genes.use = gsub("\\.2","", rownames(E13_new))
genes.unique = intersect(rownames(Count.mat), rownames(E13_new))
genes.dup = setdiff(rownames(E13_new), genes.unique)
duplicate_genes = unique(gsub("\\.[1-8]", "", genes.dup))
single_genes = setdiff(genes.unique, duplicate_genes)

CollapseCountMatrix = function(Count_data, single_genes, duplicate_genes){
  
  Count_data_new = Count_data[single_genes,]
  
  for (g in duplicate_genes){
    
    genes_to_collapse = c(g, grep(paste0(g,"\\.[1-8]"), rownames(Count_data), value=TRUE))
    print(genes_to_collapse)
    Count_data_new = rbind(Count_data_new, colSums(Count_data[genes_to_collapse,]))
    rownames(Count_data_new)[nrow(Count_data_new)] = g
  }
  
  return(Count_data_new)
  
}

Count_mat_E13 = CollapseCountMatrix(E13_new, single_genes, duplicate_genes)
Count_mat_E14 = CollapseCountMatrix(E14_new, single_genes, duplicate_genes)

Count.mat_E13toP0 = cbind(Count_mat_E13[rownames(Count.mat),],
                      cbind(Count_mat_E14[rownames(Count.mat),],
                      Count.mat[rownames(Count.mat),]))
saveRDS(Count.mat_E13toP0, file="../../RawData/Counts_E13toP0.rds")


######################################################
############ P5 Drop-seq + 10X data #######################
#####################################################

# Load data (DropSeq P5)
Count.mat_p5dsq = readRDS("../../RawData/p5_DropSeq_reSeq_2018.rds")
# Load Data (10X)
load("../../RawData/p5rgc_counts_Dec16.Rdata")

genes.to.use = intersect(rownames(Count.mat),rownames(Count.mat_p5dsq))

colnames(Count.mat) = gsub("p5RGC","P5Cd90RGC10X",colnames(Count.mat))
colnames(Count.mat_p5dsq) = gsub("P5MouseRGC","P5Cd90RGC",colnames(Count.mat_p5dsq))

Count.mat_p510x = Count.mat
plot(rowMeans(Count.mat_p510x[genes.to.use,]),rowMeans(Count.mat_p5dsq[genes.to.use,]))

Count.mat_p5 = cbind(Count.mat_p5dsq[genes.to.use,], Count.mat_p510x[genes.to.use,] )

save(Count.mat_p5, file="../../RawData/Counts_P5.rds")
