"""

Module for annotating cell types in the AnnData object.

"""

import scanpy as sc
import pandas as pd

def find_marker_genes(adata, groupby='leiden', method='wilcoxon', verbose = True):
    """
    Find marker genes for each cluster in the AnnData object

    Parameters:
    adata (AnnData): The input AnnData object.
    groupby (str): The column in adata.obs to use for grouping the data (e.g., 'leiden').
    method (str): The method to use for finding marker genes (e.g., 'wilcoxon').

    Returns:
    DataFrame: A DataFrame containing the marker genes for each cluster.
    """

    if verbose:
        print(f"Finding marker genes for {groupby} groups...")

    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)

    if verbose:
        print("Marker gene identification completed.")

    return adata

def assign_cell_types(adata, cell_type_map, cluster_col='leiden', new_col='cell_type'):
    """
    Map cluster IDs to cell type labels
    
    Parameters
    ----------
    cell_type_map : dict mapping cluster ID strings to cell type names
    cluster_col   : column in adata.obs containing cluster IDs
    new_col       : name for the new cell type column
    """
    adata.obs[new_col] = adata.obs[cluster_col].map(cell_type_map)
    print(f"Cell type distribution:")
    print(adata.obs[new_col].value_counts())
    return adata


def score_gene_set(adata, gene_list, score_name, verbose=True):
    """
    Compute a per-cell activity score for a set of genes
    Useful for scoring pathways, cell states, or disease signatures
    
    Parameters
    ----------
    gene_list  : list of gene names to score
    score_name : name for the new score column in adata.obs
    """
    # only keep genes that actually exist in the dataset
    valid_genes = [g for g in gene_list if g in adata.var_names]

    if verbose:
        print(f"Scoring {score_name}: "
              f"{len(valid_genes)}/{len(gene_list)} genes found")

    sc.tl.score_genes(adata, valid_genes, score_name=score_name)
    return adata


def export_marker_table(adata, output_path='../results/marker_genes.csv',
                        n_genes=50):
    """
    Export top marker genes per cluster to a CSV file
    For README tables and supplementary results
    """
    result = adata.uns['rank_genes_groups']
    groups = list(result['names'].dtype.names)

    rows = []
    for group in groups:
        for i in range(min(n_genes, len(result['names'][group]))):
            rows.append({
                'cluster':       group,
                'gene':          result['names'][group][i],
                'score':         result['scores'][group][i],
                'pval_adj':      result['pvals_adj'][group][i],
                'logfoldchange': result['logfoldchanges'][group][i]
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"Marker table saved to {output_path}")
    return df