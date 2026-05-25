"""
Module for clustering analysis of single-cell RNA-seq data using Scanpy. This module provides functions to run PCA, KNN, UMAP, and Leiden clustering on an AnnData object, as well as a function to run the full clustering pipeline.

"""

import scanpy as sc

def run_pca(adata, n_comps=50, verbose = True):
    """
    Run PCA on the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    n_comps (int): Number of principal components to compute.

    """
    if verbose:
        print("Running PCA...")

    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)

    if verbose: 
        print("PCA completed.")

    return adata

def run_KNN(adata, n_neighbors=10, n_pcs=50, verbose = True):
    """
    Run KNN on the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    n_neighbors (int): Number of neighbors to use for KNN graph construction.
    n_pcs (int): Number of principal components to use for KNN graph construction.

    """
    if verbose:
        print("Running KNN...")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    if verbose:
        print("KNN completed.")
    return adata

def run_UMAP(adata, n_neighbors=10, min_dist=0.3, verbose = True):
    """
    Run UMAP on the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    n_neighbors (int): Number of neighbors to use for UMAP graph construction.
    min_dist (float): Minimum distance parameter for UMAP.

    """
    if verbose:
        print("Running UMAP...")

    sc.tl.umap(adata, n_neighbors=n_neighbors, min_dist=min_dist)

    if verbose:
        print("UMAP completed.")

    return adata

def run_leiden(adata, resolution=0.5, verbose = True):
    """
    Run Leiden clustering on the AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object.
    resolution (float): Resolution parameter for Leiden clustering.

    """
    if verbose:
        print("Running Leiden clustering...")

    sc.tl.leiden(adata, resolution=resolution)

    if verbose:
        print("Leiden clustering completed.")

    return adata

def run_clustering_pipeline(adata, n_pcs=50, n_neighbors=10, min_dist=0.3, resolution=0.5, verbose = True):
    """
    Run the full clustering pipeline: PCA, KNN, UMAP, Leiden clustering.

    Parameters:
    adata (AnnData): The input AnnData object.
    n_pcs (int): Number of principal components to compute for PCA and use for KNN.
    n_neighbors (int): Number of neighbors to use for KNN and UMAP graph construction.
    min_dist (float): Minimum distance parameter for UMAP.
    resolution (float): Resolution parameter for Leiden clustering.

    """
    adata = run_pca(adata, n_comps=n_pcs, verbose=verbose)
    adata = run_KNN(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, verbose=verbose)
    adata = run_UMAP(adata, n_neighbors=n_neighbors, min_dist=min_dist, verbose=verbose)
    adata = run_leiden(adata, resolution=resolution, verbose=verbose)

    return adata