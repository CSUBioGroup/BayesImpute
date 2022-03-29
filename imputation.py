
# coding: utf-8

# In[ ]:



# coding: utf-8

# In[ ]:

import numpy as np

def find_cluster_cell_idx(l, label):
    '''
        find cells whose label == l
        return bool_idx
    '''
    return label==l 

 
def find_cluster_samples(idx,X):
    '''
        find samples of cells from idx
        X: n_features * n_samples
    '''
    X_samples=list(map(lambda x: np.sum(x != -1),X[:,idx]))
    return X_samples

  
def find_cluster_gene_mean(idx,X):
    '''
        find mean, std of cells from idx
        X: n_features * n_samples
    '''
    X_bar=list(map(lambda x:(x[x>=0]).mean(),X[:,idx]))
    return X_bar


def find_othercluster_cell_idx(l,label):
    '''
        find cells whose label != l
        return bool_idx
    '''
    return label!=l


def find_othercluster_gene_mean_std(idx, X):
    '''
        find mean, std of cells from idx
        X: n_features * n_samples
    '''
    theta=list(map(lambda x:(x[x>=0]).mean(),X[:,idx]))
    tau=list(map(lambda x:(x[x>=0]).std(),X[:,idx]))
    return theta,tau



def imputation_for_cell(cell_idx, X, labels, global_std, cluster_samples, same_cluster_stats, other_cluster_stats):
    '''
        X: n_features * n_samples
    '''
    global_std = global_std
    l = labels[cell_idx]
    same_mean = same_cluster_stats[l]
    other_mean = other_cluster_stats[l][0] 
    other_std  = other_cluster_stats[l][1] 
    gene_idx = np.arange(X.shape[0])

    
    def imputation_for_gene(gene_idx):
        # print(gene_idx, type(global_std), len(global_std), global_std[gene_idx])
        subcell_num = cluster_samples[l][gene_idx]
        delta_2 = global_std[gene_idx]**2
        
        tau_2 = other_std[gene_idx]**2
   
        denom = delta_2 + subcell_num*tau_2
   
        v = ((subcell_num*tau_2)/ denom * same_mean[gene_idx])+ (delta_2 / denom * other_mean[gene_idx])
  
        return v
    
    return list(map(imputation_for_gene, gene_idx))
