# scRNA-seq Immune Cell Profiling Pipeline

End-to-end single-cell RNA sequencing analysis pipeline characterizing 
8 immune cell populations from 2,700 PBMCs using Python and Scanpy.

## Pipeline
- Notebook 1: Quality control and filtering
- Notebook 2: Normalization and feature selection
- Notebook 3: Dimensionality reduction and clustering
- Notebook 4: Marker gene detection and cell type annotation

## Results
8 distinct immune populations identified including CD4/CD8 T cells, 
B cells, monocytes, NK cells, dendritic cells, and megakaryocytes.

## Requirements
- Python 3.10
- Scanpy, AnnData, umap-learn, leidenalg, pandas, NumPy, Matplotlib, Seaborn

## Setup
conda create -n scrna python=3.10
conda activate scrna
pip install scanpy anndata umap-learn leidenalg plotly seaborn jupyterlab