# IMPORTANT: requires python 3.8
import os
import sys
import random
import scanpy as sc
import numpy as np
import pandas as pd
import scSLAT
from scSLAT.model import run_SLAT_multi
import time
from memory_profiler import memory_usage


def load_sample(file_list):
    # function to load anndata files 
    # and add them to a list
    adata_list = []
    for file in file_list:
        print(file)
        adata_list.append(sc.read_h5ad(file))
    return adata_list


def nscluster(in_adata):
    # function to perform Leiden clustering
    # using input anndata file
    
    # Copying counts from main layer to "counts" layer
    in_adata.layers["counts"] = in_adata.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(in_adata)
    # Logarithmize the data
    sc.pp.log1p(in_adata)
    # PCA
    sc.tl.pca(in_adata)
    # Nearest neighbor graph construction
    sc.pp.neighbors(in_adata)
    # UMAP
    sc.tl.umap(in_adata)
    # Leiden clustering
    sc.tl.leiden(in_adata, n_iterations = -1)
    # copying normalised counts from main layer to "logcounts" layer
    in_adata.layers["logcounts"] = in_adata.X.copy()
    # replacing main layer data with "counts" layer
    in_adata.X = in_adata.layers["counts"].copy()
    return in_adata


def slat_run(file_list, cos = 0.3):
    # Function to run slat alignment
    # The first slice is treated as reference
    # and clustering is performed on it.
    # Remaining slices are aligned to the 
    # reference and cluster labels are borrowed 
    # from it.
    
    adata_list = load_sample(file_list)
    adata_list[0] = nscluster(adata_list[0])
    #### multislice runs - Aligning multiple datsets
    matching_list, zip_res = run_SLAT_multi(adata_list, cos_cutoff = cos, n_jobs = 1) 
    #### borrowing labels
    adata_list[0].obs['clust_ref'] = adata_list[0].obs['leiden'].copy()
    
    for ad_num in range(len(adata_list)-1):
        ref_ad = adata_list[ad_num]
        query_ad = adata_list[ad_num + 1]
        match_index = matching_list[ad_num][1]
        clust_ref = ref_ad.obs['clust_ref'].copy()
        if len(match_index) < query_ad.shape[0]:
            print("Not all matches found for query", ad_num + 1)
            for i in range(query_ad.shape[0]):
                if i not in matching_list[ad_num][0]:
                    match_index = np.insert(match_index, i, 1000000)
            clust_query = []
            for j in range(len(match_index)):
                if match_index[j] == 1000000:
                    clust_query.append("NA")
                else:
                    clust_query.append(clust_ref[ref_ad.obs_names[match_index[j]]])
            clust_query = pd.Series(clust_query)
            clust_query.index = query_ad.obs_names.copy()
        else:
            clust_query = clust_ref[ref_ad.obs_names[match_index[0:]]]
            clust_query.index = query_ad.obs_names.copy()
        query_ad.obs['clust_ref'] = clust_query
    
    return adata_list


def slat_mem_wrapper(file_list, cos = 0.3):
    start_time = time.time()
    mem_usage, adata_list = memory_usage((slat_run, (file_list, cos)), retval=True, max_usage=True)
    end_time = time.time()

    print(f"Maximum memory usage: {mem_usage} MiB")
    print(f"Total running time: {end_time - start_time} seconds")
    adata_list.append({"mem":mem_usage, 
                       "tic":start_time,
                       "toc":end_time})
    return adata_list
