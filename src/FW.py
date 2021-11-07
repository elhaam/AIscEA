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


      
def process_cisTopic_on_orig(rep_cisTopic, f_original, out_dir, save=False, transpose=True, str_to_int=False):    
    # 1) Adding the clusters of cisTpoic rep3 to the clusters of acutal rep3 anndata
    # original h5ad file
    rep = sc.read(out_dir + f_original)
    if str_to_int:
        rep.var.reset_index(drop=True, inplace=True)
        rep.obs.reset_index(drop=True, inplace=True)
    if transpose == True:
        rep = rep.T
    cisTopic_clusters = pd.DataFrame(rep_cisTopic.obs["leiden"])
    rep.obs['leiden'] = list(cisTopic_clusters["leiden"])
    rep.obs['leiden'] = rep.obs['leiden'].astype('category')
    #cisTopic_clusters['leiden'].astype('category')

    # 2) Project
    rep.obsm['X_tsne'] = rep_cisTopic.obsm['X_tsne']
    rep.obsm['X_umap'] = rep_cisTopic.obsm['X_umap']
    
    ######## Marker genes #########
    marker_genes, rep_markers_added = get_marker_genes(rep)
    return marker_genes, rep_markers_added



def get_marker_genes(adata):
    ################################### Marker genes ################################
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df_rank = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'scores', 'pvals_adj', 'logfoldchanges']})
    abs_log2fc = False
    for i in range(len(df_rank.columns)//4):
        if abs_log2fc == False:
            df_temp = df_rank.iloc[:,i*4:(i*4)+4]
            df_temp = df_temp.sort_values(by=[df_rank.columns[(i*4)+3]], ascending=False)
        else:
            df_temp = df_rank.iloc[:,i*4:(i*4)+4]
            df_temp = df_temp.iloc[(-np.abs(df_temp[df_rank.columns[(i*4)+3]].values)).argsort()]
            
        # Round logfoldchanges to 2 decimals
        df_temp[df_temp.columns[3]] = df_temp[df_temp.columns[3]].apply(lambda x: '{:.2f}'.format(round(x, 2)))
        
        # Only keep genes with p-vals threshold
        df_temp_pval = df_temp[df_temp[df_temp.columns[2]] <= 0.01]
        if i == 0:
            df_ranks = df_temp_pval.reset_index(drop=True)
        else:
            df_ranks = pd.concat([df_ranks, df_temp_pval.reset_index(drop=True)], axis=1)   
    pval_cols = [i * 4 + 2 for i in range(df_ranks.shape[1]//4)]
    df_ranks = df_ranks.drop(df_ranks.columns[pval_cols], axis=1)
    scores_cols = [i * 3 + 1 for i in range(df_ranks.shape[1]//3)]
    df_ranks_small = df_ranks.drop(df_ranks.columns[scores_cols], axis=1)
    return df_ranks, adata 


 
def match_clusters3(markers_rna, markers_ATAC, method="rank_rna_shared", show=True, normal=False, norm_m='l1', verbose =True, threshold=1, top_100=False, disp_res=False):
    '''
    Matching the clusters between RNA-seq and ATAC-seq data using various methods
    '''
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
                marker_list = df_rna.index.isin(df3_greater_1.index)
                indecies_in_rna = np.where(marker_list)[0]
                if top_100:
                    indecies_in_rna = np.asarray([i for i in indecies_in_rna if i <= 100])
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
                
            elif method == 'count':
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = len(vector_rna)
                
                n_shared_list = []
                n = 1000
                for i in range(n):
                    # Random sample of size: Number of cells in RNA having logFC above the threshold
                    sample_rna = df_rna.sample(len_rna_g_thr)
                    # Random sample of size: Number of cells in ATAC having logFC above the threshold
                    sample_atac = df_atac.sample(len_atac_g_thr)
#                     print("P-value: ", len(set(sample_rna.index) & set(sample_atac.index)))
                    n_shared_list.append(len(set(sample_rna.index) & set(sample_atac.index)))
                # Find how many times elements in n_share_list are >= res using count for that entry
                p_val.loc[str(col)+'_rna', str(col2)+'_atac'] = sum(i > res.loc[str(col)+'_rna', str(col2)+'_atac'] for i in n_shared_list) / n
            
            elif method == 'proportion':
                # Intersect of genes above the threshold
                union_markers = set(df_atac.index) & set(df_rna.index)
                print("ATAC and RNA size: ", len(df_atac.index), len(df_rna.index))
                print("Prop: ", len(vector_rna), len(union_markers))
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = len(vector_rna) / len(union_markers)
                
            elif method == 'prop_intersect':
                # shared markers
                all_markers = len(df3)
                print("ATAC and RNA size: ", len(df_atac.index), len(df_rna.index))
                print("Jaccard thr: ", len(vector_rna), all_markers)
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = len(vector_rna) / all_markers
                
            elif method == 'prop_union':
                # union of markers
                union_markers = set(df_atac.index) | set(df_rna.index)
                print("ATAC and RNA size: ", len(df_atac.index), len(df_rna.index))
                print("Jaccard thr: ", len(vector_rna), len(union_markers))
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = len(vector_rna) / len(union_markers)
                
            elif method == 'jaccard_thr':
                # union of markers above the threshold
                df_atac = df_atac[df_atac[df_atac.columns[-1]] >= threshold]
                df_rna = df_rna[df_rna[df_rna.columns[-1]] >= threshold]
                union_markers = set(df_atac.index) | set(df_rna.index)
                print("ATAC and RNA size: ", len(df_atac.index), len(df_rna.index))
                print("Jaccard thr: ", len(vector_rna), len(union_markers))
                res.loc[str(col)+'_rna', str(col2)+'_atac'] = len(vector_rna) / len(union_markers)

    if disp_res == True:           
        display(res)
#     print("P-values:") 
#     display(p_val)
    ############################ Linear assignment ################################
    try:
        cost = res.to_numpy(dtype='float') * -1
        row_ind, col_ind = linear_sum_assignment(cost)
        #cost = cost*-1
#         print(col_ind)
#         print(cost[row_ind, col_ind].sum())
        return col_ind
    except:
        return 0





    


def get_cells_cluster(adata, cluster):
    '''
    Returns list of the cells in adata in the certain cluster
    '''
    # list of the cells in this cluster 
    cluster_cells = list(adata.obs[adata.obs['leiden'] == str(cluster)]['leiden'].index)
    return cluster_cells

    

def calc_log2fc(df, cell_name, gene,  verbose):
    '''
    Returns foldchange and log2foldchange for one cell and one gene using gene expression df of a cluster
    df: DataFrame that contains the cluster expression matrix or the whole expression matrix depending on 
    the compare_mode
    '''
    if '&' in cell_name:
        cells = cell_name.split('&')  
    else:
        cells = cell_name
    df2 = df.loc[:, gene]
    expr = df2.loc[cells]
    # Remove the cell to calculate the geometric mean on the rest
    df2 = df2.drop(cells)
    
    if '&' in cell_name:
        mean_expr = np.mean(expr._get_numeric_data().mean())
    else:
        mean_expr = np.mean(expr)
    mean_rest = np.mean(df2._get_numeric_data())
    
    foldchanges = (np.expm1(mean_expr) + 1e-9) / (np.expm1(mean_rest) + 1e-9)  # add small value to remove 0's
    logfoldchanges = np.log2(foldchanges)
#     if verbose == True:
#         print(cell_name, gene, mean_rest, (np.expm1(mean_expr) + 1e-9), (np.expm1(mean_rest) + 1e-9), foldchanges, logfoldchanges)
#     print(foldchanges, logfoldchanges)
    return logfoldchanges
    
        
    
def calc_log2fc_vectors(df_cluster, df_all, verbose, method='avg'):
    '''
    Returns foldchange and log2foldchange for all cells and all genes
    df_cluster: expression matrix including cells in the cluster and their marker genes
    df_all: expression matrix including all cells and genes
    '''
    if method == 'avg':
#         display(df_cluster.head())
#         print(df_cluster.mean(axis=0))
        # expression value of a gene in a cell / avg expr value of a gene in the cluster of the cell
        df_cluster = df_cluster.div(df_cluster.mean(axis=0), axis=1)
        # Fill NA values with zero (result from when devided by zero)
        df_cluster = df_cluster.fillna(0)
#         display(df_cluster.head())
#         print(df_cluster.mean(axis=0))
#         print('.............')
        return df_cluster
    
    # Values of real logFC
    else:
        lfc_cluster = dict()
    #     lfc_df = pd.DataFrame(index=df_cluster.index, columns=df_cluster.columns)
        for cell_name in df_cluster.index:
    #         fc_list  = []
            lfc_cell = dict()
            for gene in df_all.columns:
                lfc = calc_log2fc(df_all, cell_name, gene, verbose)
    #             fc_list.append(fc)
                lfc_cell[gene] = lfc

    #         fc_dict[cell_name] = fc_list
    #         lfc_df.loc[cell_name] = lfc_cell
            lfc_cluster[cell_name] = lfc_cell
        #display(lfc_df)
        return lfc_cluster


def get_union_marker_genes(df, col_ind):
    '''
    df includes marker gene names and their attributes
    '''
    union_markers = []
    for k in range(len(col_ind)):
        markers = df[str(k) + "_n"].dropna().values
        union_markers.extend(markers)
    print(len(union_markers), len(set(union_markers)))
    union_markers = set(union_markers)
    return union_markers



def get_knn(X, n_neighbors, mode="connectivity", metric="minkowski", level=1):
    """
    Parameters
    ----------
    X : array-like of shape (n_samples, n_features) or BallTree
        Sample data, in the form of a numpy array or a precomputed BallTree.
    
    n_neighbors : int
        Number of neighbors for each sample.
    
    mode : {‘connectivity’, ‘distance’}, default=’connectivity’
        Type of returned matrix: ‘connectivity’ will return the connectivity matrix with ones and zeros, 
        and ‘distance’ will return the distances between neighbors according to the given metric.    
    
    metric : str, default=’minkowski’
        The distance metric used to calculate the k-Neighbors for each sample point. 
        The default distance is ‘euclidean’ (‘minkowski’ metric with the p param equal to 2.)  
        
    Returns
    -------
    A : sparse matrix of shape (n_samples, n_samples)
        Graph where A[i, j] is assigned the weight of edge that
        connects i to j. The matrix is of CSR format.
    """
    
    assert (mode in ["connectivity", "distance"]), "Mode argument must be either 'connectivity' or 'distance'."
    # If ‘auto’, then True is used for mode=’connectivity’ and False for mode=’distance’.
    include_self = False #if mode == "distance" else True
    
    graph_knn = kneighbors_graph(X, n_neighbors, mode, metric, include_self=include_self)
    graph_knn = graph_knn.toarray()
    if level >= 2:
        knn_2 = np.matmul(graph_knn, graph_knn)
        # A + A^2
        mat = graph_knn + knn_2
        if level == 3:
            knn_3 = np.matmul(knn_2, graph_knn)
            # A + A^2 + A^3
            mat = mat + knn_3
        # binarize mat
        mat = np.where(mat < 1, 0, 1)
        # make sure that the diagonal is all zero
        np.fill_diagonal(mat, 0)
        return mat
    return graph_knn



def make_symm(G, method='or'):
    '''
    Make the matrix symmetric
    We can make the matrix symmetric in two ways: 
    1) If node A is connected to node B then node B must also be connected to node A.
    2) We only keep the connection between two nodes if A is connected to B AND B is also connected to A.
    '''
    if method == 'or':
        for row in range(len(G)):
            for col in range(len(G[row])):
                if G[row][col] == 1:
                    G[col][row] = 1                    
    if method == 'AND':            
#         for row in range(len(G)):
#             for col in range(len(G[row])):                
#                 if G[row][col] == 1 and G[col][row] == 1:
#                     pass
#                 else:
#                     G[col][row] = 0
#                     G[row][col] = 0
        G_symm = np.maximum(G, G.transpose())
    return G_symm



def get_shared_markers(markers_rna, markers_ATAC, threshold):
    shared_markers = []
    cols_float = []
    for i in range(len(markers_rna.columns) // 3):
        cols_float.append(str(i) + "_l")
        
    markers_rna[cols_float] = markers_rna[cols_float].astype(float) 
    markers_ATAC[cols_float] = markers_ATAC[cols_float].astype(float) 
    res_index = [str(s) + '_rna' for s in list(range(len(markers_rna.columns) // 3))]
    res_col = [str(s) + '_atac' for s in list(range(len(markers_rna.columns) // 3))]
    res = pd.DataFrame(index = res_index, columns = res_col)

    common_markers_dfs = dict()
    # For all clusters in scRNA
    for col in range(len(markers_rna.columns) // 3): 
        df_rna = markers_rna[[str(col) + "_n", str(col)+ "_s", str(col)+"_l"]].dropna().set_index(str(col) + "_n")
        # For all clusters in ATAC
        col2 = col
        df_atac = markers_ATAC[[str(col2) + "_n", str(col2)+ "_s", str(col2)+"_l"]].dropna().set_index(str(col2) + "_n")
        df3 = pd.merge(df_rna, df_atac, left_index=True, right_index=True)
        shared_markers.append(df3.index)
    return shared_markers
    


def solve_FW(adj_G, adj_H, C, lambd, n_iter, verbose, show_details, opt_gamma):
    ####################### Initialize P permutation  matrix #######################
    # First initialize P as an Identity matrix , and then shuffle the rows randomly
    # This assures having exactly one 1 entry at each row and each columns
    P = np.identity(len(C))
    P = shuffle(P)

    obj_values = []
    counter = 0
    while(counter < n_iter):
        ####################### J(P): conserved interactions #######################
        # A_{P(H)} = P A_{H} P.T
        adj_p_H = np.matmul(np.matmul(P, adj_H), (P.T))
        # J(P)
        J_P = np.trace(np.matmul(adj_G.T, adj_p_H)) / 2
        
        ####################### S(P): Similarity #######################
        # S(P) = tr(PC)
        S_P = np.trace(np.matmul(P, C))

        ####################### Objective function #######################
        obj = (lambd * J_P) + ((1 - lambd) * S_P)
        obj_values.append(obj)

        ####################### Balance between J(P) and S(P) #######################
        J_grad = np.matmul(np.matmul(adj_G.T, P), adj_H)
        S_grad = C.T
        grad = (lambd * J_grad) + ((1 - lambd) * S_grad)


        ########### use linear assignment to calculate P_star (P_{n+1}) ##############
        cost = grad * -1
        row_ind, col_ind = linear_sum_assignment(cost)
        if counter % 10 == 0:
            if verbose==True:
                print(counter, lambd, obj, " Col index:", col_ind)
        # print(cost[row_ind, col_ind].sum())

        P_star = np.zeros_like(cost)
        for i in range(len(row_ind)):
            P_star[row_ind[i]][col_ind[i]] = 1
        if opt_gamma != 'opt':
            # Update P basic
            gamma = 2 / (counter + 2)
        else:
#             print('Optimizing gamma')
            ############################# Updated gamma ###############################
            # coefficient for degree 1 of gamma
            B_prime1 = (1-lambd) * np.matmul((P_star-P), C)
            B_prime2 = (lambd/2) * np.matmul(adj_G, (np.matmul(P ,np.matmul(adj_H, (P_star.T-P.T))) + np.matmul((P_star-P) ,np.matmul(adj_H, P.T))))
            B_prime = -1 * np.trace(B_prime1 + B_prime2) # flip the signs according to obj function  max -> min

            # coefficient for degree 2 of gamma
            # flip the signs according to obj function max -> min
            A_prime = -1 * np.trace((lambd/2) * np.matmul(adj_G, np.matmul((P_star-P) ,np.matmul(adj_H, (P_star.T-P.T)))))

            # Find a solution in the range [0, 1]
            if A_prime != 0:
                solution = (-1/2) * (B_prime / A_prime)
            else:
                solution = 0
                
            # Check 0, 1, and solution to see which one maximizes the objective function
            obj_solutions = dict()
            for gamma in [0, 1, solution]:
                if gamma >= 0 and gamma <=1:
                    Q = P + (gamma * (P_star - P))
                    adj_Q_H = np.matmul(np.matmul(Q, adj_H), (Q.T))
                    # J(P)
                    J_Q = np.trace(np.matmul(adj_G.T, adj_Q_H)) / 2
                    ####################### S(P): Similarity #######################
                    # S(P) = tr(PC)
                    S_Q = np.trace(np.matmul(Q, C))
                    ####################### Objective function #######################
                    obj_Q = (lambd * J_Q) + ((1 - lambd) * S_Q)
                    obj_solutions[gamma] = obj_Q
            Q_star = max(obj_solutions, key=lambda k: obj_solutions[k]) 
#             print(obj_solutions, gamma)  
        ############################## End updated gamma ##############################
        P = ((1 - gamma) * P) + (gamma * P_star)
        counter = counter + 1 
    if show_details:
        plt.plot(obj_values)
        plt.show()

    ########### linear assignment to final find permutation matrix using P ##############
    cost = P * -1
    row_ind, col_ind = linear_sum_assignment(cost)
    final_P = np.zeros_like(cost)
    for i in range(len(row_ind)):
        final_P[row_ind[i]][col_ind[i]] = 1
          
    # Calculate number of conserved edges
    # A_{P(H)} = P A_{H} P.T
    adj_p_H = np.matmul(np.matmul(final_P, adj_H), (final_P.T))

    # J(P)
    J_P = np.trace(np.matmul(adj_G.T, adj_p_H)) / 2   
#     if (J_P - int(J_P)) > 0.1:
#         for i in range(len(adj_G)):
#             unk = np.matmul(adj_G.T, adj_p_H)
#             if unk[i, i] != 0:
#                 print(unk[i, i])

    if verbose==True:
        print('Number of conserved edges: ', J_P)
        print('adj_G.T \n', adj_G.T)
        print(adj_p_H)
        print(np.matmul(adj_G.T, adj_p_H))
        # Permutation using the final P matrix solved by Frank-Wolfe
        calc_adj_H = np.matmul(np.matmul(final_P, adj_G), final_P.T)
        print('First G: \n', adj_G, '\nH:\n', adj_H, '\n')
        print('Using permutaion final_P on adj_G : \n', calc_adj_H, '\n')  
    return final_P, col_ind, J_P



def get_expanded_df(df):
    # columns
    for col in df.columns:
        if '&' in col:
            col_split = col.split('&')
            for new_col in col_split:
                df[new_col] = df[col]
            df.drop(col, axis=1)
    # Rows
    for ind in df.index:
        if '&' in ind:
            ind_split = ind.split('&')
            for new_ind in ind_split:
                df.loc[new_ind] = df.loc[ind]
            df.drop(ind, axis=0)
    return df


def get_similarity(cls, clusters_logfc_list, method='simSIP'):
    '''
    Returns the similarity between the cells in one domain (either RNA or ATAC)
    '''
    sim_all = []
    sim_all_non_scaled = []
    sim_df = pd.DataFrame()
    if method == 'simSIP':
        rank_dict = dict()
        # I used expr/avg instead of logFC
        df_cls = clusters_logfc_list[cls]
        cls_np = df_cls.to_numpy()
        for i in range(len(df_cls.values)):
            cell = df_cls.index[i]
            expr = cls_np[i,:]
            rank = get_ranking(expr)
            rank_dict[cell] = rank
        
        # Calculate similarities between each pair of cells
        for cell1, rank1 in rank_dict.items():
            for cell2, rank2 in rank_dict.items():
                # sum of geometrical mean
                simSIP = np.sum(np.sqrt(1/ (rank1 * rank2)))
    
                sim_df.loc[cell1, cell2] = simSIP
    else:
        for cell1, gene_logfc1 in clusters_logfc_list[cls].items():
            for cell2, gene_logfc2 in clusters_logfc_list[cls].items():
                count_threshold = 0
                for gene in set(gene_logfc1).intersection(set(gene_logfc2)):
                    if  gene_logfc1[gene] >= my_threshold and gene_logfc2[gene] >= my_threshold:
    #                     print(cell1, cell2, gene, gene_logfc1[gene], gene_logfc2[gene])
                        count_threshold = count_threshold + 1

                sim_df.loc[cell1, cell2] = count_threshold
    #     display(sim_df.loc[0:4, 0:3])
    
        # Scale the values
    sim_all_non_scaled.append(sim_df)
    sim_all.append(sim_df/np.max(sim_df.to_numpy()))
    return sim_all, sim_all_non_scaled     


def rm_small_cluster(adata_r, adata_a, col_ind, min_cells=150):
    '''
    Removes a cluster of cells with size of less than min_cells
    '''
    # RNA
    for cls in set(adata_r.obs['leiden']):
        if(len(adata_r.obs[adata_r.obs['leiden'] == cls]) < min_cells):
            # Remove small RNA cluster
            adata_r = rm_cls(adata_r, cls)
            # Remove its correspondence in ATAC
            adata_a = rm_cls(adata_a, col_ind[int(cls)])
            print("R1: ", cls, col_ind[int(cls)])
    
    # ATAC
    for cls in set(adata_a.obs['leiden']):
        if(len(adata_a.obs[adata_a.obs['leiden'] == cls]) < min_cells):
            # Remove small RNA cluster
            adata_a = rm_cls(adata_a, cls)
            # Remove its correspondence in ATAC
            adata_r = rm_cls(adata_r, list(col_ind).index(int(cls)))
            print("R2: ", list(col_ind).index(int(cls)), cls)
    return adata_r, adata_a


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



def scRNAseq_clustering_original(data_adrs, filtering=False, resl=0.4, highly_var=False, mus=False): 
    adata = sc.read(
    data_adrs,  
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=False)  
    print(adata) 
    if mus == False:
        adata=adata.transpose()
    print(adata)
    print(len(adata.obs.index.tolist()))
    print(len(adata.var.index.tolist()))
    
    # Remove char after "." e.g. ENSMUSG00000028713.17 (added for sciCAR)
    try:
        adata.var.index = [ind.split('.')[0] for ind in adata.var.index]
    except TypeError:
        pass


    # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    adata.var_names_make_unique()

#     filtering = False
    if filtering:
        sc.pp.filter_cells(adata, min_genes=10)
        sc.pp.filter_genes(adata, min_cells=3)
        print(adata)
        print(len(adata.obs.index.tolist()))
        print(len(adata.var.index.tolist()))
        ExpressedGene = adata.var.index.tolist()

        # Filtering low quality cells
        adata = adata[adata.obs.n_genes < 25000, :]
        print(adata)

    # Normalization and scaling
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
    sc.pp.log1p(adata)
    #sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.log1p(adata)
    #We will freeze the current state of the AnnData object, i.e. the logarithmized raw gene expression, in the a raw attribute. 
    adata.raw = adata

    if highly_var:
        # Identify highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        #sc.pp.highly_variable_genes(adata,  max_mean=10000)
    
        adata = adata[:, adata.var.highly_variable]
    print(adata)
    adata.raw = adata

    ##################################################################################
    # After here, it's the same as clustering fuction above
    # Computing neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=0, metric='cosine')
    # Embedding the neighborhood graph
    sc.tl.leiden(adata, resl)
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    

    # Figure
    rcParams['figure.figsize'] = 9, 9
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data', legend_fontsize=12, legend_fontoutline=2,frameon=False, title='UMAP', palette='gist_rainbow')#
    sc.pl.tsne(adata, color=['leiden'], legend_loc='on data', legend_fontsize=12, legend_fontoutline=2,frameon=False, title='tSNE', palette='gist_rainbow')

    # Group cells based on their cluster membership
    E = adata.to_df().reset_index()
    df = pd.DataFrame(adata.obs['leiden']).reset_index().rename(columns={'index': 'cell'})
    clusters = []
    for i in range (df['leiden'].nunique()):
        cluster = df[df['leiden'].cat.codes==i]['cell'].tolist()
        clusters.append(cluster)
        #print(f"Size of cluster {i}: {len(cluster)}")

    Es = []
    for i in range (df['leiden'].nunique()):
        Ei = E.loc[E['index'].isin(clusters[i])].set_index('index')
        print(i, ":", len(Ei.index.values.tolist()))
        Es.append(Ei)
        
    ################################### Marker genes ################################
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    #sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df_rank = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'scores', 'pvals_adj', 'logfoldchanges']})
    
    abs_log2fc = False
    
    for i in range(len(df_rank.columns)//4):
        if abs_log2fc == False:
            df_temp = df_rank.iloc[:,i*4:(i*4)+4]
            df_temp = df_temp.sort_values(by=[df_rank.columns[(i*4)+3]], ascending=False)
        else:
            df_temp = df_rank.iloc[:,i*4:(i*4)+4]
            df_temp = df_temp.iloc[(-np.abs(df_temp[df_rank.columns[(i*4)+3]].values)).argsort()]
            
#             df_temp.abs().sort_values(by=[df_rank.columns[(i*4)+3]], ascending=False)
            
#             df.iloc[(-np.abs(df['b'].values)).argsort()]
        # Round logfoldchanges to 2 decimals
        df_temp[df_temp.columns[3]] = df_temp[df_temp.columns[3]].apply(lambda x: '{:.2f}'.format(round(x, 2)))
        
        # Only keep genes with p-vals threshold
        df_temp_pval = df_temp[df_temp[df_temp.columns[2]] <= 0.01]
        if i == 0:
            df_ranks = df_temp_pval.reset_index(drop=True)
        else:
            df_ranks = pd.concat([df_ranks, df_temp_pval.reset_index(drop=True)], axis=1)
        
    
    pval_cols = [i * 4 + 2 for i in range(df_ranks.shape[1]//4)]
    df_ranks = df_ranks.drop(df_ranks.columns[pval_cols], axis=1)
    scores_cols = [i * 3 + 1 for i in range(df_ranks.shape[1]//3)]
    df_ranks_small = df_ranks.drop(df_ranks.columns[scores_cols], axis=1)
    return df_ranks, adata

