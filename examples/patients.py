import sys
# Second append the folder path in name-file.py
sys.path.insert(0, '/home/ejafari/alignment/downstream/notebooks/FW_src/')
from utils import *
from collections import Counter


def find_multi_cells(adata_k, n_clusters):
    
    X = adata_k.to_df().to_numpy()
    kmeans = KMeans(n_clusters, random_state=123).fit(X)
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
    print("After combining multi-cells: ", cluster_df.shape) 
    return cluster_df



def run_fw(X, Y, similarity,   rna_multi_clusters, cluster, n_iter, lambd, ignore=False, show_details=False, verbose=False, gamma='opt'):
    final_P, col_ind_cells, j_p = solve_FW(Y, X, similarity.values, lambd, n_iter, verbose, show_details, gamma)
    ############################# Calculate FOCCTTM score #############################
    fracs1, fracs2, fo = calc_foscttm(rna_multi_clusters[cluster], atac_multi_cells_clusters[col_ind[cluster]], col_ind_cells, similarity, cluster, ignore, show_details)
    if show_details:
        print('-----------------------------------------------------------------------')
    return final_P, col_ind_cells, fo, j_p, fracs1, fracs2



def calc_foscttm(data1, data2, col_ind_cells, similarity, cluster_num, ignore_no_cluster_match=False, show_details=False):
    fracs1 = []
    fracs2 = []
    data1_align_dict = dict()

    
    cells_in_cluster_rna = get_cells_cluster(rna, cluster_num)
    cells_in_cluster_atac = get_cells_cluster(atac_cis_on_org, col_ind[cluster_num])
    
    print('Inside calc_foscttm: ')
#     display(atac_cis_on_org.to_df().head())
    print(len(set(cells_in_cluster_rna) & set(cells_in_cluster_atac)))
    print("Check: ", len(cells_in_cluster_rna), len(cells_in_cluster_atac), col_ind[cluster_num])
    
    adata_rna_k = rna[rna.obs.loc[cells_in_cluster_rna].index].to_df() # ALL genes. TODO: explore when using ONLY MARKER genes 
    adata_atac_k = atac_cis_on_org[atac_cis_on_org.obs.loc[cells_in_cluster_atac].index].to_df() # ALL genes. TODO: explore when using ONLY MARKER genes
    print('Inside calc_foscttm, adata_atac_k: ')
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
    print(len(adata_rna_k), len(adata_atac_k), len(fracs1), frac1_mean, len(fracs2), frac2_mean, fosccttm)
    return fracs1, fracs2, fosccttm


# Remove any cluster of RNA or ATAC having less than 150 cells
def rm_small_cluster2(adata_r, adata_a, markers_r, markers_a, min_cells=150, rm_correspondence=True):
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


# Remove any cluster of RNA or ATAC having less than 150 cells
def rm_small_cluster(adata_r, adata_a, min_cells=150):
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


# Remove any cluster of RNA or ATAC having less than 150 cells
def rm_small_cluster2(adata_r, adata_a, markers_r, markers_a, min_cells=150, rm_correspondence=True):
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


# Remove any cluster of RNA or ATAC having less than 150 cells
def rm_small_cluster(adata_r, adata_a, min_cells=150):
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


def match_clusters3(markers_rna, markers_ATAC, method="dot", show=True, normal=False, norm_m='l1', verbose =True, threshold=1, top_100=False, disp_res=False):
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
                
            elif method == 'count' or method == 'rank_rna_shared':
                res_count.loc[str(col)+'_rna', str(col2)+'_atac'] = len(vector_rna)
                
                n_shared_list = []
                n = 200
                for iterr in range(n):
                    # Random sample of size: Number of cells in RNA having logFC above the threshold
                    sample_rna = df_rna.sample(len_rna_g_thr)
                    # Random sample of size: Number of cells in ATAC having logFC above the threshold
                    sample_atac = df_atac.sample(len_atac_g_thr)
#                     print("P-value: ", len(set(sample_rna.index) & set(sample_atac.index)))
                    n_shared_list.append(len(set(sample_rna.index) & set(sample_atac.index)))
                # Find how many times elements in n_share_list are >= res using count for that entry
                p_val.loc[str(col)+'_rna', str(col2)+'_atac'] = sum(i > res_count.loc[str(col)+'_rna', str(col2)+'_atac'] for i in n_shared_list) / n
            
            
                
 
        
        
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
#         display(res)
        pass
    print("P-values:") 
#     display(p_val)


    ################################## Linear assignment ##############################################
    from scipy.optimize import linear_sum_assignment
    
    try:
        cost = res.to_numpy(dtype='float') * -1
        row_ind, col_ind = linear_sum_assignment(cost)
        
        #cost = cost*-1
#         print(col_ind)
#         print(cost[row_ind, col_ind].sum())
        return col_ind, res
    except:
        return 0
    
    
def match_clusters3(markers_rna, markers_ATAC, method="dot", show=True, normal=False, norm_m='l1', verbose =True, threshold=1, top_100=False, disp_res=False):
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
                
            elif method == 'count' or method == 'rank_rna_shared':
                res_count.loc[str(col)+'_rna', str(col2)+'_atac'] = len(vector_rna)
                
                n_shared_list = []
                n = 200
                for iterr in range(n):
                    # Random sample of size: Number of cells in RNA having logFC above the threshold
                    sample_rna = df_rna.sample(len_rna_g_thr)
                    # Random sample of size: Number of cells in ATAC having logFC above the threshold
                    sample_atac = df_atac.sample(len_atac_g_thr)
#                     print("P-value: ", len(set(sample_rna.index) & set(sample_atac.index)))
                    n_shared_list.append(len(set(sample_rna.index) & set(sample_atac.index)))
                # Find how many times elements in n_share_list are >= res using count for that entry
                p_val.loc[str(col)+'_rna', str(col2)+'_atac'] = sum(i > res_count.loc[str(col)+'_rna', str(col2)+'_atac'] for i in n_shared_list) / n
            
            
                
 
        
        
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
                


#     if disp_res == True:           
#         display(res)
#     print("P-values:") 
#     display(p_val)


    ################################## Linear assignment ##############################################
    from scipy.optimize import linear_sum_assignment
    
    try:
        cost = res.to_numpy(dtype='float') * -1
        row_ind, col_ind = linear_sum_assignment(cost)
        
        #cost = cost*-1
#         print(col_ind)
#         print(cost[row_ind, col_ind].sum())
        return col_ind, res
    except:
        return 0
    
    
############################ High_var = True ######################################
input_dir = "/home/ejafari/alignment/downstream/data/Patients/mine/rna.csv"
# scRNAseq_adrs = input_dir + 'scRNAseq.csv'
# Find marker genes of scRNAseq data
markers_rna, rna = scRNAseq_clustering_original(input_dir, filtering=True, resl=1, highly_var=True, mus=False)



# Find clusters of low dimensional cisTopic result for scATACseq

input_dir = "/home/ejafari/alignment/downstream/data/Patients/mine/cisTopic/"
atac = clustering(input_dir + 'atac_topics.tsv', transpose=True, resl=0.3) # 0.5 -> 7



############################ scATAC ######################################
# 1) Put cisTopic clusters and embedding values on the original data and 2) find marker genes and their logFC
input_dir = "/home/ejafari/alignment/downstream/data/Patients/mine/"
f_original = 'pred_matrix_closest_genes_1k_downstream_strand_first_prom.csv'
markers_atac, atac_cis_on_org = process_cisTopic_on_orig(atac, f_original, input_dir, save=False, transpose=True)




# To match name of the cells in RNA and ATAC data
atac_cis_on_org.obs.index = [i.replace('.','-') for i in atac_cis_on_org.obs.index] 


import time
start = time.time()


rna, atac_cis_on_org,  markers_rna, markers_atac = rm_small_cluster2(rna, atac_cis_on_org,markers_rna, markers_atac, min_cells=200, rm_correspondence=False)



# Checkpoint

pd.set_option('display.max_columns', 30)
n, m = len(set(rna.obs['leiden'])), len(set(atac_cis_on_org.obs['leiden']))
df_aggr = pd.DataFrame(np.zeros((n,m)), columns=np.arange(0,m), index=np.arange(0,n))
for my_threshold in np.arange(0.0, 1.0, 0.1):
#     print(my_threshold)
    col_ind, res = match_clusters3(markers_rna, markers_atac, method="rank_rna_shared", verbose=False, threshold=my_threshold, top_100=True, disp_res=True)
    print(my_threshold.round(1), col_ind)
    print('$$$-------------------------------------------------------------------$$$')
    for i in range(len(col_ind)):
        df_aggr.iat[i, col_ind[i]] = df_aggr.iat[i, col_ind[i]] + 1
        
# display(df_aggr)

# Linear Assignment to get final col_ind (cluster alignment)
cost = df_aggr.to_numpy(dtype='float') * -1
row_ind, col_ind = linear_sum_assignment(cost)
print(col_ind)




# Step 1
# Values of the above dataframe which are > 5 
# Which means they have been assigned together at least half of ten times


# SRR matrix for threshold 0
my_threshold = 0
col_ind_0, SRR = match_clusters3(markers_rna, markers_atac, method="rank_rna_shared", verbose=False, threshold=my_threshold, top_100=True, disp_res=True)
print(my_threshold, col_ind_0)
print("--------------------%%--------------------")

# Check values of col_ind to see if they are > 5
for row in range(len(col_ind)):
    
    if df_aggr.loc[row, col_ind[row]] < 6:
        print(row, col_ind[row], df_aggr.loc[row, col_ind[row]], "   Removed!")
        # Remove RNA cluster and corresponding ATAC cluster 
        # for example remove BOTH RNA cluster 2 that corresponds to  ATAC cluster 1 with value of 5.0
        # RNA
        rna = rm_cls(rna, row)
#         ATAC
        atac_cis_on_org = rm_cls(atac_cis_on_org, col_ind[row])
        
    elif df_aggr.loc[row, col_ind[row]] == 6:
        # If the value is 6 it means that the mapping is not certain
        # We need to check and see how many other clusters are bigger in that row in df_aggr matrix when thr = 0
        # Get the (sum of the absolute ditance to this cluster alignment) /  (the value for this alignment in df_aggr)

        
        # SRR value of mapped clusters in the matrix when threshold = 0
        print('$$$------------------------- Value 6 -------------------------------$$$')
        SRR_mapping = SRR.iloc[row, col_ind[row]]
        print(row, col_ind[row], SRR_mapping)
        
        # Check the row of SRR matrix for RNA cluster to get all other similarity score greater than SRR_mapping
        SRR_row = list(SRR.iloc[row,:])
        print(SRR_row)
        # Subtract the original SRR for mapped cluster from all the values
        SRR_row_diff = [diff for diff in SRR_row - SRR_mapping if diff > 0] 
        print("Diff", SRR_row_diff)
        SRR_diff_sum = np.sum(SRR_row_diff)
        print("Diff sum: ", SRR_diff_sum, "--- Diff sum / SRR mapping: ", SRR_diff_sum / SRR_mapping)
        print("-----------------------------------------------------------------")
        
        if SRR_diff_sum / SRR_mapping > 0.5:
            # Remove
            print(row, col_ind[row], df_aggr.loc[row, col_ind[row]], "   Removed greater than 6!")
            # RNA
            rna = rm_cls(rna, row)
    #         ATAC
            atac_cis_on_org = rm_cls(atac_cis_on_org, col_ind[row])
            
        else:
            # Keep
            print(row, col_ind[row], df_aggr.loc[row, col_ind[row]], "   Kept although greater than 6!")
    else:
        print(row, col_ind[row], df_aggr.loc[row, col_ind[row]], "   Kept!")
        
        
# Double check the true mappings using cell labels

df_cells_common = pd.DataFrame()

for cls_rna in sorted(set(rna.obs['leiden'])):
    for cls_atac in sorted(set(atac_cis_on_org.obs['leiden'])):
#         cells_atac  = atac_cis_on_org.obs[atac_cis_on_org.obs['leiden'] == str(cls_atac)].index
        cells_atac = [i.replace('.','-') for i in atac_cis_on_org.obs[atac_cis_on_org.obs['leiden'] == str(cls_atac)].index]
        cells_rna = rna.obs[rna.obs['leiden'] == str(cls_rna)].index
        df_cells_common.loc[('r_' + str(cls_rna)), ('a_' + str(cls_atac))] = len(set(cells_atac) & set(cells_rna))
        
df_cells_common

# 


# From the above we know
col_ind = {0:3, 1:2, 2:1, 3:0}

# Taken from Darwin
import time
from collections import Counter
from sklearn.cluster import KMeans

time1 = time.time()
# Hyperparameter 1
for my_threshold in [0]:
    print(my_threshold)
#     col_ind = match_clusters4(markers_rna, markers_atac, method="count", verbose=False, threshold=my_threshold)
    
    union_markers_rna = get_union_marker_genes(markers_rna, col_ind)
    union_mrakers_atac = get_union_marker_genes(markers_atac, col_ind)
    intersect_marker_genes = (union_markers_rna & union_mrakers_atac)
    print("Intersect: ", len(intersect_marker_genes))

    # Save cells of rna and atac clusters in two seperate dictionaries
    clusters_r = dict()
    clusters_a = dict()
    atac_k_all_genes = []
    rna_k_all_genes = []

    atac_multi_cells_clusters = dict()
    rna_multi_cells_clusters = dict()  
    
   ##################################  conserved edges ##################################
    clusters_r_intersect = dict()
    clusters_a_intersect = dict()
    clusters_r_before_multi = dict()
    clusters_a_before_multi = dict()
    intersect_marker_r_a = get_shared_markers(markers_rna, markers_atac, my_threshold)
    
    # clusters of scRNA
    for rna_k, atac_k in col_ind.items():
        # aligned cluster of scATACseq
#         atac_k = col_ind[rna_k]
        print(rna_k, atac_k)   

        # marker genes in this cluster of rna
        markers_rna_k = markers_rna[str(rna_k) + "_n"].dropna().values
        cells_in_cluster_rna = get_cells_cluster(rna, rna_k)
        # DataFrame of expression matrix for each cluster
        rna_k_all_genes.append(rna[rna.obs.loc[cells_in_cluster_rna].index, :].to_df())
        # Only marker genes
        adata_rna_k = rna[rna.obs.loc[cells_in_cluster_rna].index, markers_rna_k]
        print("RNA: ", adata_rna_k.shape)

        # marker genes in this cluster of ATAC
        markers_atac_k = markers_atac[str(atac_k) + "_n"].dropna().values
        cells_in_cluster_atac = get_cells_cluster(atac_cis_on_org, atac_k)
        # DataFrame of expression matrix for each cluster
        atac_k_all_genes.append(atac_cis_on_org[atac_cis_on_org.obs.loc[cells_in_cluster_atac].index, :].to_df())
        adata_atac_k = atac_cis_on_org[atac_cis_on_org.obs.loc[cells_in_cluster_atac].index, markers_atac_k]
        print("ATAC: ", adata_atac_k.shape)

        print("------------------------before Kmeans----------------------")
        time2 = time.time()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)

        ########################################### KMeans #################################
        if adata_rna_k.shape[0] > adata_atac_k.shape[0]:
            # Perform Kmeans on RNA data using size of ATAC data
            rna_multicell_df = find_multi_cells(adata_rna_k, n_clusters=adata_atac_k.shape[0])
            rna_multi_cells_clusters[rna_k] = rna_multicell_df # Modified 
            atac_multi_cells_clusters[col_ind[rna_k]] = adata_atac_k.to_df() # Not modified

        elif adata_rna_k.shape[0] < adata_atac_k.shape[0]:
            # Perform Kmeans on ATAC data using size of RNA data
            atac_multicell_df = find_multi_cells(adata_atac_k, n_clusters=adata_rna_k.shape[0])
            atac_multi_cells_clusters[col_ind[rna_k]] = atac_multicell_df
            rna_multi_cells_clusters[rna_k] = adata_rna_k.to_df()
        
        print("------------------------after Kmeans----------------------")
        time2 = time.time()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)

        ################# Calculate logFC for each cell based on its cluster's marker genes #################
        
        
        log2fc_rna_cluster = calc_log2fc_vectors(rna_multi_cells_clusters[rna_k], rna.to_df()[intersect_marker_genes], verbose=False)
        log2fc_atac_cluster = calc_log2fc_vectors(atac_multi_cells_clusters[col_ind[rna_k]], atac_cis_on_org.to_df()[intersect_marker_genes], verbose=False)

        clusters_r[rna_k] = log2fc_rna_cluster
        clusters_a[atac_k] = log2fc_atac_cluster
        
        # Multi cells
        log2fc_rna_cluster = calc_log2fc_vectors(rna_multi_cells_clusters[rna_k], rna.to_df()[intersect_marker_r_a[rna_k]], verbose=False)
        log2fc_atac_cluster = calc_log2fc_vectors(atac_multi_cells_clusters[col_ind[rna_k]], atac_cis_on_org.to_df()[intersect_marker_r_a[rna_k]], verbose=False)

        clusters_r_intersect[rna_k] = pd.DataFrame(log2fc_rna_cluster)
        clusters_a_intersect[atac_k] = pd.DataFrame(log2fc_atac_cluster)
        print("------------------------after multicells----------------------")
        time2 = time.time()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)

    print("------------------------before similarity----------------------")  
    time2 = time.time()
    print("A >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)
    ################################## Cell level similarities between RNA and ATAC  ##################################
    sim_all = dict()
    sim_all_non_scaled = dict()

    print_p = True

    # For each cluster
    for cls_i, cls_a in col_ind.items():
        expr_vect1_greater_dict = dict()
    #         cls_a = col_ind[cls_i]
        print(cls_i, cls_a)
        sim_df = pd.DataFrame()

        for cell1 in clusters_r[cls_i].index:
            expr_vect1 = clusters_r[cls_i].loc[cell1]
            # greater than 1 -> sort
            expr_vect1_greater = expr_vect1[expr_vect1 >= 1].sort_values(ascending=False)
    #         print(cell1, expr_vect1_greater)
            expr_vect1_greater_dict[cell1] = expr_vect1_greater
            if print_p:
                print('RNA')
#                 print(expr_vect1_greater)
                print_p = False       
        time2 = time.time()
        print("A >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)

        print_p = True
        for cell2 in clusters_a[cls_a].index:
            expr_vect2 = clusters_a[cls_a].loc[cell2]

            # Find genes of ATAC that have avg expr >= 1 which are shared and >=1  in both RNA and ATAC
            expr_vect2_greater = expr_vect2[expr_vect2 >= 1]
            for cell1, expr_vect1_greater in expr_vect1_greater_dict.items():

                # Shared genes >= 1
                shared_genes_g_1 = set(expr_vect1_greater.index).intersection(expr_vect2_greater.index)

                genes_to_rank = expr_vect2_greater.loc[shared_genes_g_1].index

                 # Find rank RNA shared similairty
                marker_list = expr_vect1_greater.index.isin(expr_vect2_greater.index)
                indecies_in_rna = np.where(marker_list)[0]
                if True: #Top 100
                    indecies_in_rna = np.asarray([i for i in indecies_in_rna if i <= 100])
                if print_p:
                    print('ATAC')
#                     print(expr_vect2_greater)
                    print("#Genes shared and greater than 1 in RNA and ATAC: ", len(shared_genes_g_1))
                    print_p = False
                    print(indecies_in_rna)


                sim_df.loc[cell1, cell2] = np.sum(np.sqrt(1 / (indecies_in_rna + 1)))
        sim_all[cls_i] = sim_df/np.max(sim_df.to_numpy()) # scale

    print("------------------------after RNA-ATAC similarity----------------------")   
    time2 = time.time()
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)



    print("------------------------before similarity RNA-RNA and ATAC-ATAC ----------------------")  
    time2 = time.time()
    print("A >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)   
    sim_rna_all = dict()
    sim_atac_all = dict()
    # For each cluster    
    for cls_i, cls_a in col_ind.items():
        ################################## Cell level similarities in RNA  #######################
        sim_all1, sim_all_non_scaled11 = get_similarity(cls_i, clusters_r)
        sim_rna_all[cls_i] = sim_all1
        #sim_rna_all_non_scaled.append(sim_all_non_scaled1)
        ################################## Cell level similarities in ATAC  #######################
        sim_all1, sim_all_non_scaled1 = get_similarity(col_ind[cls_i], clusters_a)
        sim_atac_all[col_ind[cls_i]] = sim_all1
        #sim_atac_all_non_scaled.append(sim_all_non_scaled1)   

    print("------------------------before FW----------------------")   
    time2 = time.time()
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", time2 - time1)        

    #3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 45, 50, 60, 70, 80, 90
    for n_neigh in [3]:
    #         if n_neigh < 30:
    #             levs = [1, 2, 3]
    #         elif n_neigh < 50:
    #             levs = [1, 2]
    #         else:
    #             levs = [1]
        levs = [2]
        for lev in levs:

            adj_rna = dict()
            adj_atac = dict()

            # Hyperparameter 2   
            lambds = [0.8] # 0.2, 0.8
            ind = 0
            fosccttm_df = pd.DataFrame(columns=['Cluster', 'Lambda', 'fracs1_sum', 'len_fracs1', 'fracs2_sum', 'len_fracs2'])
            n_runs = 1
            #             true_shared_E_list = []
            #             fracs_sum = pd.DataFrame(columns=['cluster', 'correspondence', 'sum_fracs', 'len'])

            fracs1_all = []
            fracs2_all = []
            
            for i, atac_i in col_ind.items():
            #                 atac_i = col_ind[i]
                ind, col = clusters_r[i].index, clusters_r[i].index
                vals = get_knn(clusters_r[i], n_neigh, level=lev)
                # Symmetric KNN
                vals = make_symm(vals, method='AND')
                adj_rna[i] = pd.DataFrame(vals, index=ind, columns=col)
                ind, col = clusters_a[atac_i].index, clusters_a[atac_i].index
                vals = get_knn(clusters_a[atac_i], n_neigh, level=lev)
                # Symmetric KNN
                vals = make_symm(vals, method='AND')
                adj_atac[atac_i] = pd.DataFrame(vals, index=ind, columns=col)



            #             for i in range(len(adj_rna)): 
                X = adj_rna[i].copy()
                Y = adj_atac[atac_i].copy()
                print(X.shape, Y.shape)

                X_mult_sim = np.multiply(X, sim_rna_all[i][0])
                Y_mult_sim = np.multiply(Y, sim_atac_all[atac_i][0])


                print("************************************ Cluster ", i, "************************************ ")
                for l in lambds:

                    for n_iter in range(n_runs):
                        final_P, col_ind_cells, frac, j_p, fracs_list1, fracs_list2 = run_fw(X_mult_sim.values, Y_mult_sim.values, sim_all[i], rna_multi_cells_clusters, i, 40, l) # gamma='opt' default
                        fracs1_all = fracs1_all + fracs_list1
                        fracs2_all = fracs2_all + fracs_list2


                        fosccttm_df.loc[0 if pd.isnull(fosccttm_df.index.max()) else fosccttm_df.index.max() + 1] = [int(i), l, np.sum(fracs_list1), len(fracs_list1), np.sum(fracs_list2), len(fracs_list2)]
            #                         fracs_sum.loc[0 if pd.isnull(fracs_sum.index.max()) else fracs_sum.index.max() + 1] = [i, col_ind[i], frac, len_fo]

            #                 np.save(str(i) + "_" + str(col_ind[i]) + str(n_neigh) + 'knn_marker_genes_thr_' + str(my_threshold) + '_level_' + str(lev) + "_" +  str(n_iter) +  '_fracs1.npy', fracs_list1) # save
            #                 np.save(str(i) + "_" + str(col_ind[i])  + str(n_neigh) + 'knn_marker_genes_thr_' + str(my_threshold) + '_level_' + str(lev) +  "_" + str(n_iter) +   '_fracs2.npy', fracs_list2) # save 
            #display(fosccttm_df)
            print(n_neigh, my_threshold)
            #             fosccttm_df.to_csv(str(n_neigh) + 'knn_marker_genes_thr_' + str(my_threshold) + '_level_' + str(lev) + '.csv', index=True, header=True, sep='\t')



# X on Y
print((fosccttm_df['fracs1_sum'].sum() / (fosccttm_df['len_fracs1'].sum() - 1)) / fosccttm_df['len_fracs1'].sum())
# fosccttm_df['fracs1_sum'].sum() / fosccttm_df['len_fracs1'].sum()

# Y on X
print((fosccttm_df['fracs2_sum'].sum() / (fosccttm_df['len_fracs2'].sum() - 1)) / fosccttm_df['len_fracs2'].sum())
# fosccttm_df['fracs2_sum'].sum() / fosccttm_df['len_fracs2'].sum()


# Save fracs_all1 and fracs_all2

print(len(fracs1_all), len(fracs2_all))

import pickle
out_dir = "/home/ejafari/alignment/downstream/notebooks/Patients/FW_res/"

with open(out_dir + 'fracs1.pickle', 'wb') as handle:
    pickle.dump(fracs1_all, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
with open(out_dir + 'fracs2.pickle', 'wb') as handle:
    pickle.dump(fracs2_all, handle, protocol=pickle.HIGHEST_PROTOCOL)
            
end = time.time()
duration = end - start
print("Duration:", duration)
