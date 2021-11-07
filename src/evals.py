# Last updated: Aug, 2021
__author__ = "Elham Jafari"
__email__ = "ejafari@indiana.edu"   

import umap
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
import collections
from collections import namedtuple
import numpy as np
import pandas as pd
import scanpy as sc
import umap
# import umap.plot
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.preprocessing import StandardScaler, normalize



def get_foscttm(rna, atac):
    """ 
    This function gets the embedded scRNAseq and the embedded scATACseq and returns 
        the FOSCTTM (fraction of samples closer than the true match) values.
  
    Parameters: 
    rna (2D array): embedded scRNAseq data 
    atac (2D array): embedded scATACseq data 
    
    Returns: 
    foscttm (float): fraction of samples closer than the true match
  
    """
    
    foscttm = []
    if len(rna) != len(atac):
        print('Wrong input, size mismatch.')
        return False
    
    RNA_ATAC = namedtuple("dist", ["RNA_ind", "ATAC_ind"])
    # Compute the distance between all scRNAseq embedded points and all scATACseq embedded points
    for i in range(len(rna)):
        dist_dict = {}
        for j in range(len(atac)):
            # Save distances between a single RNA point to all ATAC points in a dictionary
            dist_dict[str(j)] = dist(rna[i], atac[j])

            
        # Sort the distance dict by distance between RNA point vs all ATAC points
        #print('Sorting the distances to get FOSCTTM:')
        dist_sorted = sorted(dist_dict.items(), key=lambda x: x[1])

        # Find the number of samples closer than the true match for each cell
        location = 0
        for tup in dist_sorted:
            if tup[0] == str(i):
                #print(tup, location)
                foscttm.append(location / len(rna))
                #print(foscttm)
            else:
                location = location + 1
                
    plt.scatter(range(len(rna)), sorted(foscttm), marker='.', s=4) 
    plt.xlabel("Cells (sorted)")
    plt.ylabel("FOSCTTM (The lower, the better)")
    plt.show()
    
    # Fraction of samples closer than the true match
    return dist_sorted, foscttm


def dist(a, b):
    """
    Returns Euclidean ditance
    """
    return np.sqrt(np.sum((a - b) ** 2, axis=0))


def get_n_closer(cell1, cell2, adata_k, cells_in_cluster, show_details):
    if cell1 in cells_in_cluster:
        if cell1 == cell2:
            n_closer = 0
            return n_closer
        else:    
            true_match = adata_k.loc[cell1, :].values
            cell2_vector = adata_k.loc[cell2, :].values
            data_without = adata_k.drop([cell1, cell2])
            mat = data_without.values
            euc_dist = np.sqrt(np.sum(np.square(np.subtract(mat, cell2_vector)), axis=1))
            dist_true_fw = np.sqrt(np.sum(np.square(np.subtract(true_match, cell2_vector))))
            n_closer = sum(i < dist_true_fw for i in euc_dist)
#             print(n_closer)
            return n_closer
    else:
        if show_details:
            print("Details: ", cell1)
        return -1