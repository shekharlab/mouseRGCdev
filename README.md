# mouseRGCdev

A summary of the main scripts used in this paper

1. `DataPrep/0_RGCdev_consolidate_counts.R` - Concatenate count matrices from different ages. Saves `RawData/Counts_E13toP0.rds` and `Counts_P5.rds`. P5 Data was kept separate as it was primarily generated using Drop-seq

2. `Clustering/1_RGCdev_combine.R` - Combined clustering analysis to define major cell classes, separate RGCs from the rest, and save S4 objects corresponding to RGCs at each time point 
