# mouseRGCdev

A summary of the main scripts used in this paper

1. `0_DataPrep.R` - This script combines count matrices from various experiments and tags each cell (column) with an ID that indicates its sample of origin. The script saves combined count matrices (genes x cells) as two separate files `RawData/Counts_E13toP0.rds` and `RawData/Counts_P5.rds`. P5 Data was kept separate as it was primarily generated using Drop-seq

2. `1_InitialClustering_Filtering.R` - Combined clustering analysis to define major cell classes, separate RGCs from the rest, and save S4 objects corresponding to RGCs at each time point. Also generates plots displayed in Figures 1, S1.

3. `2_ClarkGiudice_comparison.R` - Joint analysis of E14, E16 and P0 cells in this dataset with corresponding whole retinal scRNA-seq datasets.
     *  Clark et al., *Neuron*, 2019. Count matrices corresponding to E14, E16 and P0 were downloaded from `https://github.com/gofflab/developing_mouse_retina_scRNASeq`
     * Giudice et al., *Development*, 2019. Count matrix (E15.5 cells) was obtained from the author upon request. 
