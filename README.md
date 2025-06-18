# Cell_X_Gene

This repository contains code for analyzing and visualizing single-cell gene expression data using the [cellxgene](https://cellxgene.cziscience.com/) ecosystem. The project enables exploration of transcriptomic profiles across cell types and conditions, with tools for preprocessing, normalization, and downstream visualization.

## 📦 Features

- Load and explore `.h5ad` (AnnData) formatted single-cell datasets
- Filter doublet and low-quality cells from analysis
- (Optional) integrate multiple datasets from different
- Run low-dimensional transformation of data (normalization, PCA, UMAP, clustering)
- Run cell type annotation using differential expression and gene set enrichment analysis
- Create html markdown report of data 
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

2. **Set up a virtual environment (recommended):**
```
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows
```

3. **Install dependencies:**
```
pip install -r requirements.txt
```
5. **Run the analysis**
```
python3 /home/ewalsh/Cell_X_Gene/main.py download -s TEST \
-id 856c1b98-5727-49da-bf0f-151bdb8cb056 \
-org "Homo sapiens"

python3 /home/ewalsh/Cell_X_Gene/main.py lowD -s TEST

python3 /home/ewalsh/Cell_X_Gene/main.py GSEA -s TEST \
-tissues Retina \
-species "Homo sapiens" \
-db_species "Human"

python3 /home/ewalsh/Cell_X_Gene/main.py Report -s TEST
```
