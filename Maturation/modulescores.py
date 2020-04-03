import numpy as np
import scanpy as sc
from anndata import AnnData
from typing import Union, Optional, Tuple, Collection, Sequence, Iterable

def module_score(
    adata:AnnData,
    genes_use: list,
    score_name: Optional[str] = None,
    verbose: bool = True):
    
    """\
    Compute module scores for all cells in adata as described in methods of RGC-dev paper.
    
    
    Parameters
    ----------
    adata
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    genes_use
        list of genes in module of interest
    score_name
        Name endowed to the module score to be computed
        e.g. "Mod1"
    verbose
        Inform user of fraction of module genes that are in adata
        
    Returns
    -------
    adata with a new .obs called score_name
    
    """
    
    if (score_name==None):
        score_name = str(input("Provide a name for this score (no spaces): "))
        
    genes_use0 = genes_use
    genes_use = list(set(genes_use).intersection(adata.var_names))#genes that are both in module and `adata`
    
     
    if (len(genes_use) == 0):
        raise ValueError("Error : Must provide a list of genes that are present in the data")
        
    
    if (verbose):
        if(len(genes_use0) > len(genes_use)):
            n = len(genes_use0) - len(genes_use)
            print("Note that", n, "of the", len(genes_use0), "genes in your module do not exist in the data set." )
    
    
    
    adata_score = adata.copy()
    adata_score = adata[:,genes_use]
    
    counts_modgenes = adata_score.X.toarray() #all cells, module genes
    counts_all = adata.X.toarray() #all cells, all genes
    scores = np.mean(counts_modgenes, axis=1) - np.mean(counts_all, axis=1) #(row means of counts_modgenes ) - (row means of counts_all)
    
    adata.obs[score_name] = scores
    
    return genes_use


#Test 1 for operation, not accuracy
h5ad_path = '../forFLE/FLE.h5ad'
FLE = sc.read_h5ad(h5ad_path)

with open('mod1.txt', 'r') as f:
    mod1 = [line.strip() for line in f]

x = module_score(FLE, mod1, "Mod1Score")


print(FLE.obs)
print(len(x))



print("~~~~~~~~~~`Space between tests~~~~~~~~~~")

#Test 2 for accuracy

mod1_indices = []
for i in range(0,len(mod1)):
    
    if (mod1[i] in FLE.var_names):
        mod1_indices.append(FLE.var_names.get_loc(mod1[i]))

counts_modgenes = np.zeros((101)) #all cells, module genes
counts_all = FLE.X.toarray() #all cells, all genes

for i in range(0,len(mod1_indices)):
    j =mod1_indices[i]
    counts_modgenes[i] = counts_all[0,j]
a = round(np.mean(counts_modgenes) - np.mean(counts_all[0,:]),3) 
b = round(FLE.obs['Mod1Score'][0],3)
c = round(a/b,3)

if (c == 1):
    print("Test passed with ", "a =", a, "and b =", b)
    
else:
    print("Test failed")
