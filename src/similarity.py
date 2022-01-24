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

from collections import Counter

import sys
sys.path.insert(0, '/home/ejafari/alignment/Git/src/')
from utils import *
from rmCls import *
from AIscEA import *
from evals import *



def match_clusters(markers_rna, markers_ATAC, method="dot", show=True, normal=False, norm_m='l1', verbose =True, threshold=1, top_100=False, disp_res=False):
#     print("Here", method)

    cols_float = []
    for i in range(len(markers_rna.columns) // 3):
        cols_float.append(str(i) + "_l")
    markers_rna[cols_float] = markers_rna[cols_float].astype(float) 
    cols_float = []
    for i in range(len(markers_ATAC.columns) // 3):
        cols_float.append(str(i) + "_l")
    markers_ATAC[cols_float] = markers_ATAC[cols_float].astype(float) 

    res_index = [str(s) + '_rna' for s in list(range(len(markers_rna.columns) // 3))]
    res_col = [str(s) + '_atac' for s in list(range(len(markers_ATAC.columns) // 3))]
    res = pd.DataFrame(index = res_index, columns = res_col)
    res_count = pd.DataFrame(index = res_index, columns = res_col)
    p_val = pd.DataFrame(index = res_index, columns = res_col)

    common_markers_dfs = dict()

    # For all clusters in scRNA
    for col in range(len(markers_rna.columns) // 3): 
        df_rna = markers_rna[[str(col) + "_n", str(col)+ "_s", str(col)+"_l"]].dropna().set_index(str(col) + "_n")
        # For all clusters in ATAC
        for col2 in range(len(markers_ATAC.columns) // 3):
            df_atac = markers_ATAC[[str(col2) + "_n", str(col2)+ "_s", str(col2)+"_l"]].dropna().set_index(str(col2) + "_n")
            # Number of cells in RNA data having logFC above the threshold
            len_rna_g_thr = (df_rna[df_rna.columns[-1]] >= threshold).sum()
            len_atac_g_thr = (df_atac[df_atac.columns[-1]] >= threshold).sum()
            if verbose:
                print("Number of marker genes with values greater than thr in ATAC and RNA: " , len_atac_g_thr, len_rna_g_thr)
            df3 = pd.merge(df_rna, df_atac, left_index=True, right_index=True)
            
            if verbose == True:
                print(col, col2, " shared marker genes no thr:", len(df3), "ATAC & RNA markers \wo thr:", df_atac.shape, df_rna.shape)
#                 display(df3)

                
                
                plt.figure(figsize=(18,8))
                plt.plot(np.array(df3[df3.columns[1]]), '.', markersize=8, color="red", alpha=0.6)
                plt.plot(np.array(df3[df3.columns[3]]), '.', markersize=8, color="blue", alpha=0.6)
                plt.ylabel('LogFC')
                plt.title('RNA and ATAC : ' + str(col) + ', ' + str(col2))
                plt.axhline(y=threshold, color='grey', linestyle='-')
                
                plt.show()
                
            # Find the marker genes that are greater than one in both RNA and ATAC cluster pair that we are cpmparing
            df3_greater_1 = df3[(df3[df3.columns[-1]] >= threshold) & (df3[df3.columns[1]] >= threshold)]
            

            
            
            if verbose == True:
                print("df3_greater_than_threshold size: ", df3_greater_1.shape)
            # sum up the logFC values of the RNA in the df3_greater_1
            vector_rna = df3_greater_1[df3_greater_1.columns[1]]
            vector_atac = df3_greater_1[df3_greater_1.columns[3]]
            

            if method == 'simSIP':
#                 print(method)
                rank_rna = get_ranking(df3_greater_1[df3_greater_1.columns[1]])
                rank_atac = get_ranking(df3_greater_1[df3_greater_1.columns[3]])
#                 print("RNA: ", rank_rna)
#                 print(rank_atac)
                # sum of geometrical mean
                simSIP = np.sum(np.sqrt(1/ (rank_rna * rank_atac)))
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = simSIP
            
            if method == 'rank_rna_shared':
#                 display(df3_greater_1)
#                 display(df_rna)

                marker_list = df_rna.index.isin(df3_greater_1.index)
                indecies_in_rna = np.where(marker_list)[0]
#                 print(type(indecies_in_rna))
                if top_100:
                    indecies_in_rna = np.asarray([i for i in indecies_in_rna if i <= 100])
#                 print(col, col2)
#                 print(indecies_in_rna)
#                 print(type(indecies_in_rna))
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = np.sum(np.sqrt(1 / (indecies_in_rna + 1)))

                
            if len(df3_greater_1) != 0 and normal == True:
                 try:
                    vector_rna = normalize(vector_rna.to_numpy().reshape(1, -1), norm=norm_m, axis=1)
                    vector_atac = normalize(vector_atac.to_numpy().reshape(1, -1), norm=norm_m, axis=1)
#                     print(vector_rna, vector_atac)
                 except ValueError:
                    # If size of the vectors is 0
                    res.loc[str(col)+'_rna', str(col2)+'_atac'] = 0
            
            if len(df3_greater_1) == 0:
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = 0
            elif method == "dot":
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = np.dot(vector_rna, vector_atac.T)
            elif method == 'cosine':
                try:
                    res.loc[str(col)+'_rna', str(col2)+'_atac'] = cosine_similarity(vector_rna.to_numpy().reshape(1, -1), vector_atac.to_numpy().reshape(1, -1))
                except AttributeError:
                    res.loc[str(col)+'_rna', str(col2)+'_atac'] = cosine_similarity(vector_rna.reshape(1, -1), vector_atac.reshape(1, -1))
                
            elif method == 'corr':
#                 print(np.corrcoef(vector_rna.to_numpy().reshape(1, -1), vector_atac.to_numpy().reshape(1, -1)))
                if len(vector_rna) < 2:
                    res.loc[str(col)+'_rna', str(col2)+'_atac'] = 0
                else:
                    res.loc[str(col)+'_rna', str(col2)+'_atac'] = np.corrcoef(vector_rna.to_numpy().reshape(1, -1), vector_atac.to_numpy().reshape(1, -1))[0][1]
            


    if disp_res == True:           
        display(res)


    ################################## Linear assignment ##############################################
    from scipy.optimize import linear_sum_assignment
    
    try:
        
        # Remove all_zero row/columns
        # Remove all-zero columns
        res = res[res.columns[res.sum(axis=0) != 0]]
        # Remove all-zero rows
        res = res[(res.sum(axis=1) != 0)] 
        
        
        cost = res.to_numpy(dtype='float') * -1
        row_ind, col_ind = linear_sum_assignment(cost)
        return row_ind, col_ind, res
        

    except:
        return row_ind, col_ind, res
    
    
def find_multi_cells(adata_k, n_clusters):
    
    X = adata_k.to_df().to_numpy()
    kmeans = KMeans(n_clusters, random_state=123, max_iter=20).fit(X)
    ################################## find multi-cell clusters #################################
    # find number of occurance for each cluster label
    cluster_df = adata_k.to_df()
    to_be_dropped = []
    counter_dict = Counter(kmeans.labels_)
    for i_label in range(len(kmeans.labels_)):
        
        if counter_dict[i_label] > 1:
#           print(i_label, counter_dict[i_label])

            # Add the multi-cells with their center values in KMeans
            indices = [i for i, x in enumerate(kmeans.labels_) if x == i_label]
            inds = cluster_df.index[indices]
            to_be_dropped.extend(inds)
            
            # new index that contains the name of all cells in the multi-cell
            new_index = '&'.join(inds)
#             print(new_index)
            # Mean value for the multi-cell is exactly the KMeans center
            cluster_df.loc[new_index] = kmeans.cluster_centers_[i_label]
    # drop the cells with frequencies greater than 1
#     print(to_be_dropped)
    cluster_df = cluster_df.drop(to_be_dropped)
#     print("After combining multi-cells: ", cluster_df.shape) 
    return cluster_df



def cluster_mapping(markers_rna, markers_atac, rna, atac_cis_on_org):    
    '''
    Ruturns cluster mapping after checking linear assignment for 10 different threshold of log2fc of marker genes
    -------------------
    Input:
        markers_rna (dataframe): marker genes of RNA
        markers_atac (dataframe): marker genes of ATAC
        rna (AnnData): RNA data
        atac_cis_on_org (AnnData): atac_cis_on_org data
    Output:
        col_ind (list of integers): Cluster mapping
    
    '''
    
    n, m = len(set(rna.obs['leiden'])), len(set(atac_cis_on_org.obs['leiden']))
    df_aggr = pd.DataFrame(np.zeros((n,m)), columns=np.arange(0,m), index=np.arange(0,n))
    for my_threshold in np.arange(0.0, 1, 0.1):
#         print(my_threshold)

        row_ind, col_ind, res = match_clusters(markers_rna, markers_atac, method="rank_rna_shared", verbose=False, threshold=my_threshold, top_100=True, disp_res=False)
#         display(res)
#         print(my_threshold.round(1), col_ind)
        



        
        # Calculate p-values
        iter_i = 1000
        ratio_greater = dict()
        for i in range(iter_i):
            perm_SRR = []
            cluster_map_S = dict()
            
            #permute values in col_ind to get random selections
            permutation = np.random.permutation(list(col_ind))
        #     print(permutation)
            # Find the values in SRR0 associated with permutation vector
            for ind in range(len(permutation)):
                perm_SRR.append(res.iloc[ind, permutation[ind]])
                # Save cluster mapping in a dict in a form {tuple (RNA_cls, ATAC_cls): S value}
                RNA_cls = int(res.index[ind].split("_")[0])
                ATAC_cls = int(res.columns[col_ind[ind]].split("_")[0])
                cluster_map_S[RNA_cls, ATAC_cls] = res.iloc[ind, col_ind[ind]]
        #     print(perm_SRR)

            # Sort the SRR velues in the  permutation vector
            perm_SRR_sorted = np.sort(perm_SRR)
            # Sort the dict by values of S
            cluster_map_S_sorted = dict(sorted(cluster_map_S.items(), key=lambda item: item[1]))
#             print(perm_SRR_sorted)

        #     for all entries of col_ind: compare i-th element with i-th element of perm_SRR_sorted
            for ind in range(len(perm_SRR_sorted)):
        #         print(SRR0.iloc[ind, col_ind[ind]], perm_SRR_sorted[ind])
                if list(cluster_map_S_sorted.values())[ind] >= perm_SRR_sorted[ind]:
        #             print("T")
                    try:
                        ratio_greater[list(cluster_map_S_sorted.keys())[ind]] += 1
                    except KeyError:
                        ratio_greater[list(cluster_map_S_sorted.keys())[ind]] = 1


#         print(col_ind)
#         print(ratio_greater)
        ratio_greater = {k: v / iter_i for k, v in ratio_greater.items()}
        p_val_new = {k: 1 - v for k, v in ratio_greater.items()}
#         a = {k: v / total for k, v in a.iteritems()}
#         print("P-value: ", p_val_new)

#         print("$$$$$$$ \n ")
        
#         col_ind_cor = dict()
#         for i in sorted(set(rna.obs['leiden'])):
#             # List of corresponding clusters
#             col_ind_cor[int(i)] = col_ind[int(i)] 
#         print(col_ind_cor)
#         print("$$$$$$$ \n ")

        for row in range(len(col_ind)):
            RNA_cls = int(res.index[row].split("_")[0])
            ATAC_cls = int(res.columns[col_ind[row]].split("_")[0])
            if p_val_new[RNA_cls, ATAC_cls] <= 0.01:
                # keep
#                 print(RNA_cls, ATAC_cls, df_aggr.loc[RNA_cls, ATAC_cls], "   Kept statistically significant!")
#                 for i in range(len(col_ind)):
                df_aggr.iat[RNA_cls, ATAC_cls] = df_aggr.iat[RNA_cls, ATAC_cls] + 1
            else:
                pass
                # Remove
#                 print(RNA_cls, ATAC_cls, df_aggr.loc[RNA_cls, ATAC_cls], "   Removed  p-value > 0.01!")
                # RNA
#                 rna = rm_cls(rna, row)
        #         ATAC
#                 atac_cis_on_org = rm_cls(atac_cis_on_org, col_ind[row])     
        
#     try:            
#         display(df_aggr)
#     except:
#         pass
    return df_aggr, rna, atac_cis_on_org


def extract_mapped_clusters(rna, markers_rna, atac_cis_on_org, markers_atac, p_val_count=6):
    '''
    After mapping the clusters together and calculating p-values: keep only cluster mapping with good p-values
    '''
    df_aggr, rna, atac_cis_on_org = cluster_mapping(markers_rna, markers_atac, rna, atac_cis_on_org)
    # Keep clusters with siginificant p-values (threshold = 0.01) for at least 6 log2FC out of 10
    df = df_aggr[df_aggr >= p_val_count]
    try:
        col_ind_cor = {k:v for (k, v) in list(df.stack().index)}
#         print("\nMapped clusters between Domain1 and Domain2: ", col_ind_cor)

#         print("\nRemoving non-significant clusters.")
        # Remove clusters with bad p-value which are not in col_ind_corr
        for cls_rna in sorted(set(rna.obs['leiden'])):

            if int(cls_rna) not in col_ind_cor.keys():
                # Remove
                rna = rm_cls(rna, cls_rna)

        for cls_atac in sorted(set(atac_cis_on_org.obs['leiden'])):
            if int(cls_atac) not in col_ind_cor.values():
                # Remove
                atac_cis_on_org = rm_cls(atac_cis_on_org, cls_atac)
        return col_ind_cor, rna, atac_cis_on_org
    # When no mapped clusters overall
    except IndexError:
        return None, None, None
