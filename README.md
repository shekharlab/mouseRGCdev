# mouseRGCdev

A summary of the main scripts used in the paper "Diversification of multipotential postmitotic mouse retinal ganglion cell precursors into discrete types" by Shekhar et al., 2021.

1. `0_DataPrep.R` - This script combines count matrices from various experiments and tags each cell (column) with an ID that indicates its sample of origin. The script saves combined count matrices (genes x cells) as two separate files `RawData/Counts_E13toP0.rds` and `RawData/Counts_P5.rds`. P5 Data was kept separate as it was primarily generated using Drop-seq

2. `1_InitialClustering_Filtering.R` - Combined clustering analysis to define major cell classes, separate RGCs from the rest, and save S4 objects corresponding to RGCs at each time point. Also generates plots displayed in Figures 1, S1.

3. `2_ClarkGiudice_comparison.R` - Joint analysis of E14, E16 and P0 cells in this dataset with corresponding whole retinal scRNA-seq datasets.
     *  Clark et al., *Neuron*, 2019. Count matrices corresponding to E14, E16 and P0 were downloaded from `https://github.com/gofflab/developing_mouse_retina_scRNASeq`
     * Giudice et al., *Development*, 2019. Count matrix (E15.5 cells) was obtained from the author upon request. 
     
4. `3_rgc[E13,E14,E16,P0,P5]_clustering.R` - Dimensionality and clustering analysis of RGCs at E13, E14, E16, P0 and P5 respectively. S4 objects corresponding to each time point are saved in each case. 

5. `4_Clusteredness_analysis/Clusteredness_analysis_Fig2.R` - Analysis of transcriptomic diversity and separatedness of RGC clusters at different time points. 


The data used in the study is also publicly available in processed form here:

https://drive.google.com/drive/folders/1zLxKKt7xLiVE1NdOvKrIPZ3OjvG418yl?usp=sharing

Data are organized by age. In each folder, there should be the following files (as an example for E13),

- `rgcE13_processed.mtx` : log-normalized expression matrix (genes x cells)
- `barcodes_E13.tsv` : Cell barcodes
- `features_E13.tsv` : Gene names (M. musculus)
- `rgcE13_metadata.tsv` : Metadata containing SampleID, Louvain Cluster, and columns for OT-calculated trajectory values for each terminal RGC type C1-C45. For example, the column "C14_traj" represents the fate associations of E13 cells with type C14. For more information on computing trajectories, see here:
https://broadinstitute.github.io/wot/tutorial/#notebook-3-inferring-long-range-temporal-couplings 

 Please feel free to email kshekhar@berkeley.edu for any questions. 
