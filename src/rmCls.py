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
from scipy.optimize import linear_sum_assignment
from sklearn.metrics.pairwise import cosine_similarity, linear_kernel     
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from sklearn.neighbors import kneighbors_graph
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
from scipy.optimize import linear_sum_assignment
from scipy.stats import gmean
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, '/home/ejafari/alignment/Git/src/')
from utils import *
from FW import *
from evals import *
from similarity import *



def rm_tiny_cluster(adata_r, adata_a, markers_r, markers_a, min_cells=150, rm_correspondence=True):
    '''
    Removes a cluster of cells with size of less than min_cells
    '''
    # RNA
    for cls in set(adata_r.obs['leiden']):
        if(len(adata_r.obs[adata_r.obs['leiden'] == cls]) < min_cells):
            # Remove small RNA cluster
            adata_r = rm_cls(adata_r, cls)
            # Remove its correspondence in ATAC
            if rm_correspondence:
                adata_a = rm_cls(adata_a, col_ind[int(cls)])
                print("R1: ", cls, col_ind[int(cls)])
            else:
                try:
                    markers_r.drop([str(cls)+'_n', str(cls)+'_s', str(cls)+'_l'], inplace=True, axis=1)
                    print("RNA:", cls)
                except:
                    print("RNA removed already:", cls)
    
    # ATAC
    for cls in set(adata_a.obs['leiden']):
        if(len(adata_a.obs[adata_a.obs['leiden'] == cls]) < min_cells):
            # Remove small RNA cluster
            adata_a = rm_cls(adata_a, cls)
            # Remove its correspondence in ATAC
            if rm_correspondence:
                cls_r = list(col_ind).index(int(cls))
                adata_r = rm_cls(adata_r, cls_r)
                print("R2: ", cls_r, cls)
            else:
                markers_a.drop([str(cls)+'_n', str(cls)+'_s', str(cls)+'_l'], inplace=True, axis=1)
                print("ATAC:", cls)
    return adata_r, adata_a, markers_r, markers_a
    

def rm_cls(adata, cls):
    '''
     Removes the  cells that belong to cluster cls
    '''
    cls = str(cls)
    # Get the list of cells in the cluster
    dropped_cells = adata.obs[adata.obs['leiden'] == cls].index
    drop_list = adata.to_df().index.isin(dropped_cells)
    # Remove the cluster that is smaller than min_cells
    adata = adata[~drop_list, :]
    return adata

