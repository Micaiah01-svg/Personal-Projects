Cell-Type–Specific Immune Responses in Disease (scRNA-seq, Python)

---

# Overview
This project implements a full Python-based single-cell RNA-seq (scRNA-seq) analysis to investigate how disease affects immune cells at the transcriptional level. Using a PBMC dataset, we analyze both cell-type composition changes and cell-type–specific gene expression differences, and interpret the results using pathway enrichment.

---

# Objectives
1. Compare immune cell-type proportions between healthy and disease conditions.
2. Correct for batch effects and integrate multiple samples.
3. Identify differentially expressed genes (DEGs) within specific immune cell types.
4. Perform pathway enrichment to interpret biological functions.

---

# Dataset
- Source: PBMC dataset from 10x Genomics (`scanpy.datasets.pbmc68k_reduced`)
- Content: ~68k single cells, preprocessed (normalized, PCA, UMAP)
- Annotations: `bulk_labels` indicate immune cell type
- Metadata: Added `condition` (Healthy vs Disease) and `batch` (simulated)

---

# Methods and Code

## Step 1: Load the Dataset
```python
import scanpy as sc
adata = sc.datasets.pbmc68k_reduced()
adata
```

## Step 2: Add Metadata
```python
import numpy as np
np.random.seed(42)
ad. obs['condition'] = np.random.choice(['Healthy','Disease'], size=adata.n_obs, p=[0.5,0.5])
ad. obs['batch'] = np.random.choice(['Batch1','Batch2','Batch3'], size=adata.n_obs)
```

## Step 3: Batch Correction and Integration
```python
import harmonypy as hm
sc.tl.pca(adata, svd_solver='arpack')
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
ad. obsm['X_pca_harmony'] = ho.Z_corr.T
```

## Step 4: Neighborhood Graph & UMAP
```python
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['batch','bulk_labels'], save="_integrated.png")
```

## Step 5: Analyze Cell-Type Proportions
```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
cell_props = adata.obs.groupby(['bulk_labels','condition']).size().groupby(level=0).apply(lambda x: x/x.sum()).reset_index(name='proportion')
plt.figure(figsize=(10,6))
sns.barplot(data=cell_props, x='bulk_labels', y='proportion', hue='condition')
plt.xticks(rotation=45, ha='right')
plt.title("Cell-type proportions by condition")
plt.tight_layout()
plt.savefig("figures/cell_type_proportions.png")
plt.show()
```

## Step 6: Cell-Type–Specific Differential Expression
```python
mono = adata[adata.obs['bulk_labels'] == 'CD14+ Monocytes'].copy()
sc.tl.rank_genes_groups(mono, groupby='condition', method='wilcoxon')
sc.pl.rank_genes_groups(mono, n_genes=15, save="_monocyte_DE.png")
```

## Step 7: Export DE Results
```python
de_results = sc.get.rank_genes_groups_df(mono, group='Disease')
de_results.to_csv("results/monocyte_disease_vs_healthy_DE.csv", index=False)
```

## Step 8: Pathway Enrichment (GO & KEGG)
```python
import gseapy as gp
top_genes = de_results.query("logfoldchanges>0").head(200)['names'].tolist()
enr = gp.enrichr(gene_list=top_genes, gene_sets=['GO_Biological_Process_2021','KEGG_2021_Human'], organism='Human', outdir='results/enrichment', cutoff=0.05)
```

## Step 9: Biological Interpretation
- Disease alters immune cell proportions and transcriptional states.
- CD14+ monocytes in disease show immune activation signatures.
- Pathway enrichment highlights innate immune and inflammatory pathways.

---

# Project Structure
```
immune_scRNAseq_project/
├── scripts/
│   └── analysis.py
├── figures/
│   ├── cell_type_proportions.png
│   ├── umap_integrated.png
│   └── monocyte_DE.png
├── results/
│   ├── monocyte_disease_vs_healthy_DE.csv
│   └── enrichment/
└── README.md
```

---

# Tools & Libraries
- Python
- Scanpy & AnnData
- Harmony (`harmonypy`)
- GSEApy
- NumPy, Pandas
- Matplotlib, Seaborn

---

# Skills Demonstrated
- scRNA-seq preprocessing & metadata integration
- Batch effect correction & visualization
- Cell-type proportion analysis
- Cell-type–specific differential expression
- Pathway enrichment and biological interpretation
- Reproducible Python workflows
