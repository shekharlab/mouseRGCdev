# mouseRGCdev

A summary of the main scripts used in this paper

1. `DataPrep/0_RGCdev_consolidate_counts.R` - This script combines count matrices from various experiments and tags each cell (column) with an ID that indicates its sample of origin. The script saves combined count matrices (genes x cells) as two separate files `RawData/Counts_E13toP0.rds` and `Counts_P5.rds`. P5 Data was kept separate as it was primarily generated using Drop-seq

2. `Clustering/1_RGCdev_combine.R` - Combined clustering analysis to define major cell classes, separate RGCs from the rest, and save S4 objects corresponding to RGCs at each time point 
