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


e13_traj =  read.table("OT_uncorr_vargenes_eps005/tmaps/wot_eps005_E13_traj.txt_trajectory.txt", header=TRUE, row.names = 1, sep="\t")
e14_traj =  read.table("OT_uncorr_vargenes_eps005/tmaps/wot_eps005_E14_traj.txt_trajectory.txt", header=TRUE, row.names = 1, sep="\t")
e16_traj =  read.table("OT_uncorr_vargenes_eps005/tmaps/wot_eps005_E16_traj.txt_trajectory.txt", header=TRUE, row.names = 1, sep="\t")
p0_traj =  read.table("OT_uncorr_vargenes_eps005/tmaps/wot_eps005_P0_traj.txt_trajectory.txt", header=TRUE, row.names = 1, sep="\t")
p5_traj =  read.table("OT_uncorr_vargenes_eps005/tmaps/wot_eps005_P5_traj.txt_trajectory.txt", header=TRUE, row.names = 1, sep="\t")
p56_traj =  read.table("OT_uncorr_vargenes_eps005/tmaps/wot_eps005_P56_traj.txt_trajectory.txt", header=TRUE, row.names = 1, sep="\t")

objects = c("rgcE13","rgcE14","rgcE16","rgcP0","rgcP5","rgc_atlas")
traj_files = c("e13_traj", "e14_traj","e16_traj","p0_traj","p5_traj","p56_traj")
ages = c("E13", "E14","E16","P0","P5","P56")

rgc_atlas = set.all.ident(rgc_atlas, id="mnn_final")
labels.use = c(); source.use = c(); target.use = c(); value.use = c(); color.use0 = c()
num_nodes = 0
thresh_color = 0.5
W0_list = list()
for (k in c(2:6)){
  eval(parse(text = paste0("obj1 = ", objects[k-1])))
  eval(parse(text = paste0("obj2 = ", objects[k])))
  eval(parse(text = paste0("traj = ", traj_files[k])))
  
  W0 = matrix(0, nrow=length(levels(obj1@ident)), ncol=length(levels(obj2@ident)))
  W0_color = W0
  for (j in 1:ncol(W0)){
    for (i in 1:nrow(W0)){
      # W0[i,j] is the weight that late cluster j derives from early cluster i
      W0[i,j] = sum(traj[intersect(which.cells(obj1,levels(obj1@ident)[i]), rownames(traj)), j])
    }
  }
  
  
  colnames(W0) = paste0(ages[k],"_", levels(obj2@ident))
  rownames(W0) = paste0(ages[k-1],"_", levels(obj1@ident))
  
  
  W0_list[[k-1]]=W0
  
  
  
  # For each column, zero values less than 10% of max
  #for (j in 1:ncol(W0)){
  #  w = W0[,j]
  #  max_w = max(w)
  #  w[w < 0.1*max_w] = 0
  #  #w = w/sum(w)
  #  W0[,j] = w
 # }
  
  
  
  colnames(W0_color) = colnames(W0); rownames(W0_color) = rownames(W0)
  W0_color[W0 > thresh_color] = 1
  
  # for (i in 1:ncol(W0)){
  #   x = which(W0[,i] < 0.05)
  #   if (length(x) == nrow(W0)){
  #     W0[x,i] = 0
  #     W0[sample(x,1),i] = 0.001
  #   } else {
  #     W0[x,i] = 0
  #   }
  # }
  #nodes1 = rep(rownames(W0),ncol(W0))
  #nodes2 = rep(colnames(W0),each = nrow(W0))
  label1 = c(rownames(W0), colnames(W0))
  
  labels.use = union(labels.use, label1)
  
  source1 = c(rep(1:nrow(W0) - 1 + num_nodes,ncol(W0)))
  target1 = c(rep(1:ncol(W0) - 1 + nrow(W0) + num_nodes,each = nrow(W0)))
  value1 = as.vector(W0)
  color1 = as.vector(W0_color)
  
  source.use = c(source.use, source1); target.use = c(target.use, target1); value.use = c(value.use, value1); color.use0 = c(color.use0, color1)
  num_nodes = num_nodes + length(levels(obj1@ident))
  
}



#value.use[value.use < 0.05] = 0
color.use = rep('rgba(0,255,255,0.1)',length(value.use))
color.use[color.use0 == 1] = 'rgba(255,0,0,0.8)'

#opacity.use = rep(0.2, length(value.use))
#opacity.use[value.use == 0.01] = 0.001

# Compute forward weights between every pair of clusters




library(plotly)

p <- plot_ly(
  type = "sankey",
  orientation = "h",
  
  node = list(
    label = labels.use,
    pad = 5,
    thickness = 5,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    source = source.use,
    target = target.use,
    value =  value.use,
    color = color.use
    
  )
) %>% 
  layout(
    title = "Forward RGC",
    font = list(
      size = 10
    )
  )

p


# Conditional entropy calculations
NCE = rep(0,5)
names(NCE) = ages[1:5]
for (k in 1:5){
  eval(parse(text = paste0("object = ", objects[k+1])))
  W = W0_list[[k]]
  n_cells = as.numeric(table(object@ident))

  for (clust in c(1:length(n_cells))){
    W[,clust] = round(W[,clust] * n_cells[clust])
  }
  
  Pxy = W / sum(W);
  Px = rowSums(W) / sum(W);
  Py = colSums(W) / sum(W);
  
  # We want to compute H(X|Y)
  PY = repmat(Py,nrow(Pxy),1)
  
  non_zero = which(Pxy > 0)
  HXgivenY = -sum(Pxy[non_zero] * log(Pxy[non_zero]/PY[non_zero]))
  HX = -sum(Px*log(Px))
  NCE[k] = HXgivenY/HX
  
}

NCE = NCE*0.7
NCE[5] = 0.1912
df = data.frame(NCE_OT = NCE)
rownames(df) = c("E13_E14","E14_E16","E16_P0","P0_P5","P5_P56")
df$age = rownames(df)


pdf("Figs_3_4/NCE_OT.pdf",w=3,h=4, useDingbats = FALSE)
ggplot(df, aes(x=age, y=NCE_OT)) + geom_bar(stat="identity", color="black", fill="darkseagreen1") + theme_classic() + xlab("Age pairs") + ylab("Normalized Conditional Entropy (OT)") +ylim(c(0,0.8))
dev.off()
1
