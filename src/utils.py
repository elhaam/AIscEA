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
from AIscEA import *
from evals import *
from rmCls import *
from similarity import *



def plot_dim_reduction(data1, data2):
    print('Performing tSNE:')
    embed1 = TSNE(n_components=2, random_state=0).fit_transform(data1)
    embed2 = TSNE(n_components=2, random_state=0).fit_transform(data2)
    myplot(embed1, embed2)
    print('\n')

    print('Performing UMAP:')
    embed1 = umap.UMAP().fit_transform(data1)
    embed2 = umap.UMAP().fit_transform(data2)
    myplot(embed1, embed2)
    print('\n')
    
    
      

def myplot(data1, data2, name=''):   
    plt.subplot(121)
    plt.title('scATACseq')
    plt.scatter(data1[:,0], data1[:,1], c="red", alpha=0.4, s=5)

    plt.subplot(122)
    plt.title('scRNAseq')
    plt.scatter(data2[:,0], data2[:,1], c="blue", alpha=0.4, s=5)
    plt.show()

    plt.scatter(data1[:,0], data1[:,1], c="red", alpha=0.4, s=5)
    plt.scatter(data2[:,0], data2[:,1], c="blue", alpha=0.4, s=5)
    plt.show()
    plt.savefig(name)
    
    
def atac_rna_plot(atac_file, rna_file):
    pred_matrix_atac = pd.read_csv(atac_file, delimiter='\t', index_col=0)
    pred_matrix_rna = pd.read_csv(rna_file, delimiter='\t')

    display(pred_matrix_atac.iloc[0:3, 0:5])
    display(pred_matrix_rna.iloc[0:3, 0:5])

    # convert chromosome regions in pred_matrix of rna to their associated gene names
    pred_matrix_rna.rename(gene_dict, inplace=True)


    # # Intersection of the same gene names between closest gene pred matrix and scRNAseq 
    intrsct_genes = set(pred_matrix_rna.index) & set (pred_matrix_atac.index)
    
    # How many in common gene names between closest gene pred matrix and scRNAseq 
    print("in common gene names between closest gene pred matrix and scRNAseq: ", len(intrsct_genes))
    
    # RNA - Closest gene ATAC
    print("RNA - Closest gene in ATAC: ", len(set(pred_matrix_rna.index) - set (pred_matrix_atac.index)))

    # Closest gene ATAC - RNA
    print("ATAC - Closest gene in RNA: ", len(set(pred_matrix_atac.index) - set (pred_matrix_rna.index)))
    print('--------------------------------------------------------')

    # Keep only in common genes in columns
    pred_matrix_atac = pred_matrix_atac.T[intrsct_genes]
    pred_matrix_rna = pred_matrix_rna.T[intrsct_genes]
    print(pred_matrix_rna.shape, pred_matrix_atac.shape)

    dist_sorted, foscttm = get_foscttm(pred_matrix_atac.values, pred_matrix_rna.values)
    print("Mean FOSCTTM: ", np.mean(foscttm))

    plot_dim_reduction(pred_matrix_rna.values, pred_matrix_atac.values)
    
    return pred_matrix_atac, pred_matrix_rna


def rename_gene_to_chr(df):
    gene_dict = dict()
    for i in range(df.shape[0]):
        new_name = "chr" + str(i) + ":" + str(2*i+1) + "-" + str(2*i+2)
        gene_dict[new_name] = df.index[i] 
        df.rename({df.index[i]: new_name}, inplace=True)
    return df, gene_dict



def marker_genes(adata):
    
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



def compare_atac_rna(rna, atac, markers_r, markers_a):
    print("\n--------------- Compare atac vs rna ----------------\n")
    counter = 0
    not_match = []
    for i in range(len(rna.obs['leiden'])):
    #     counter_cluster = 0
        if (rna.obs['leiden'][i] == atac.obs['leiden'][i]):
            counter = counter + 1
    #         counter_cluster = counter_cluster + 1
        else:
            not_match.append((rna.obs['leiden'][i], atac.obs['leiden'][i]))
    #     print("Number of cells that match for each cluster: \n", i, counter_cluster)
    print("-----------Number of cells with the same clusters (out of 1047 cells):------------\n", counter)
    print("\n\n")
    print(not_match)
    
#     print("\n ---------Correlations between marker genes of scATACseq clusters vs sRNAseq clusters------------ \n")
#     print('Top 50 marker genes (when in common) of scATACseq in scRNAseq data')
#     for col in range(len(markers_a.columns)//3):
#         # Top 50
#         genes = markers_a.iloc[:, col * 3].dropna().values[0:50]
#         genes = [gene for gene in genes if gene in rna.var.index]
#         print(col, genes)
#         plots_genes = ['leiden'] + list(genes)
#         sc.pl.umap(adata_rna, color=plots_genes, color_map=matplotlib.cm.Greens)
       
    print('--------------------------------------------------------------------\nTop 50 marker genes of scRNAseq in scATACseq data')
    for col in range(len(markers_r.columns)//3):
        # Top 50
        genes = markers_r.iloc[:, col * 3].dropna().values[0:50]
        genes = [gene for gene in genes if gene in atac.var.index]
        print(col, genes)
        plots_genes = ['leiden'] + list(genes)
        sc.pl.umap(adata_ATAC, color=plots_genes, color_map=matplotlib.cm.Greens)
        
        
        
def unit_normalize(data, norm="l2", sample_wise=True):
    '''
    Default norm used is l2-norm. More options: "l1" and "max"
    sample_wise: If we independently normalize each sample or we independently normalize each feature
    From SCOT: https://github.com/rsinghlab/SCOT
    '''
    assert (norm in ["l1","l2","max"]), "Norm argument has to be either one of 'max', 'l1', or 'l2'."
    if sample_wise==True:
        axis=1
    else:
        axis=0
    return normalize(data, norm=norm, axis=axis) 



