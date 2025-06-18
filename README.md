# Cell_X_Gene

This repository contains code for analyzing and visualizing single-cell gene expression data using the [cellxgene](https://cellxgene.cziscience.com/) ecosystem. The project enables exploration of transcriptomic profiles across cell types and conditions, with tools for preprocessing, normalization, and downstream visualization.

## 📦 Features

- Load and explore `.h5ad` (AnnData) formatted single-cell datasets
- Subset by cell type, condition, or gene expression
- Generate UMAP and PCA plots for visualizing clusters
- Plot gene expression distributions across clusters
- Export filtered datasets and summary statistics

## 🧪 Technologies Used

- Python 3.x
- `scanpy` for single-cell analysis
- `matplotlib` and `seaborn` for visualization
- `pandas` and `numpy` for data manipulation

## 🚀 Getting Started

1. **Clone the repo:**

   ```bash
   git clone https://github.com/ewalsh3963/Cell_X_Gene.git
   cd Cell_X_Gene
