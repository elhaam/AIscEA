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
import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.insert(0, '/home/ejafari/alignment/Git/src/')
from utils import *
from AIscEA import *
from evals import *
from rmCls import *
from similarity import *



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
    
    











def calc_foscttm(rna, atac_cis_on_org, data1, data2, col_ind_cells, similarity, cluster_num, col_ind, ignore_no_cluster_match=False, show_details=False):
    fracs1 = []
    fracs2 = []
    data1_align_dict = dict()

    
    cells_in_cluster_rna = get_cells_cluster(rna, cluster_num)
    cells_in_cluster_atac = get_cells_cluster(atac_cis_on_org, col_ind[cluster_num])
    
#     print('Inside calc_foscttm: ')
#     display(atac_cis_on_org.to_df().head())
#     print(len(set(cells_in_cluster_rna) & set(cells_in_cluster_atac)))
#     print("Check: ", len(cells_in_cluster_rna), len(cells_in_cluster_atac), col_ind[cluster_num])
    
    adata_rna_k = rna[rna.obs.loc[cells_in_cluster_rna].index].to_df() # ALL genes. TODO: explore when using ONLY MARKER genes 
    adata_atac_k = atac_cis_on_org[atac_cis_on_org.obs.loc[cells_in_cluster_atac].index].to_df() # ALL genes. TODO: explore when using ONLY MARKER genes
#     print('Inside calc_foscttm, adata_atac_k: ')
#     display(adata_atac_k.head())
    
    for i in range(len(col_ind_cells)):
        # If cell name is multi-cell (contains &) then split it 
        cell1_list = similarity.index[i].split('&')  # RNA cell
        cell2_list = similarity.columns[col_ind_cells[i]].split('&') # ATAC cell that's matched to RNA cell check prev code
        data1_align_dict[tuple(cell1_list)] = cell2_list # cell2_list is a list of size 1 if not a multi-cell
    
    
                
    for cell1_list, cell2_list in data1_align_dict.items():
        
        
        for cell1 in cell1_list:
            cell1_fracs = []
            for cell2 in cell2_list:
                ############################## Find distance in data2 (atac) ##############################
                # Cluster level distance
                n_closer = get_n_closer(cell1, cell2, adata_atac_k, cells_in_cluster_atac, show_details)
                if n_closer < 0 and ignore_no_cluster_match == False:
                    # Search the whole global space (global distance)
                    # If the cell does NOT have a true alignment in the aligned cluster
                    #Calculate distance of the node to all other cells GLOBALLLY
                    n_closer = get_n_closer(cell1, cell2, atac_cis_on_org.to_df(), atac_cis_on_org.to_df().index, show_details)
                
                # This makes sure not to include a cell if it doesn't exist in both domains
                if n_closer >= 0:
                    cell1_fracs.append(n_closer)
                
            # Mean fracs when multi-cells if correspondence exists
            if len(cell1_fracs) > 0:
                fracs1.append(np.mean(cell1_fracs))
#                 print("cell1_fracs: ", cell1_fracs)
#                 print(fracs1)

        # Find distance in data 1 (RNA)
        for cell2 in cell2_list:
            cell2_fracs = []
            for cell1 in cell1_list:
                ############################## Find distance in data1 (rna) ##############################
                # Cluster level distance
                n_closer = get_n_closer(cell2, cell1, adata_rna_k, cells_in_cluster_rna, show_details)
                if n_closer < 0 and ignore_no_cluster_match == False:
                    # Search the whole global space (global distance)
                    # If the cell does NOT have a true alignment in the aligned cluster
                    #Calculate distance of the node to all other cells GLOBALLLY
                    n_closer = get_n_closer(cell2, cell1, rna.to_df(), rna.to_df().index, show_details)
                if n_closer >= 0:
                    cell2_fracs.append(n_closer)
            # Mean fracs when multi-cells if correspondence exists
            if len(cell2_fracs) > 0:
                fracs2.append(np.mean(cell2_fracs))



    frac1_mean = np.mean(fracs1) / (len(rna) - 1)
    frac2_mean = np.mean(fracs2) / (len(atac_cis_on_org) - 1)
    
    fosccttm = np.mean([frac1_mean, frac2_mean])
#     print(len(adata_rna_k), len(adata_atac_k), len(fracs1), frac1_mean, len(fracs2), frac2_mean, fosccttm)
    return fracs1, fracs2, fosccttm, data1_align_dict


def get_ranking(arr):
    '''
    Returns rankings in an array with possible duplicate values
    '''
    u, v = np.unique(arr, return_inverse=True)
    ranks = (np.cumsum(np.bincount(v)) - 1)[v]
#     if len(arr) < 6:
#         print(ranks)
    max_rank = np.max(ranks) + 1
    ranks = [max_rank-i for i in ranks]
    return np.array(ranks)


def get_cells_cluster(adata, cluster):
    '''
    Returns list of the cells in adata in the certain cluster
    '''
    # list of the cells in this cluster 
    cluster_cells = list(adata.obs[adata.obs['leiden'] == str(cluster)]['leiden'].index)
    return cluster_cells

