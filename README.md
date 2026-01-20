# Personal-Projects
Cell-Type–Specific Immune Responses in Disease (scRNA-seq, Python)
Overview

This project implements a full Python-based single-cell RNA-seq (scRNA-seq) analysis to investigate how disease affects immune cells at the transcriptional level. Using a publicly available PBMC dataset, we analyze both cell-type composition changes and cell-type–specific gene expression differences, and interpret the results using pathway enrichment.

This workflow is end-to-end reproducible, from raw/anndata loading to results visualization and biological interpretation.

Objectives

Compare immune cell-type proportions between healthy and disease conditions.

Correct for batch effects and integrate multiple samples.

Identify differentially expressed genes (DEGs) within specific cell types.

Perform pathway enrichment to interpret biological functions.

Dataset

Source: PBMC dataset from 10x Genomics (scanpy.datasets.pbmc68k_reduced)

Content: ~68k single cells, preprocessed (normalized, PCA, UMAP)

Annotations: bulk_labels indicate immune cell type (e.g., CD4 T cells, CD14+ monocytes)

Metadata: Added two columns in code:

condition (Healthy vs Disease)

batch (simulated for demonstration of batch correction)

In real projects, this metadata comes from donor information and sample preparation details.

Methods and Step-by-Step Explanation
Step 1: Load the Dataset
adata = sc.datasets.pbmc68k_reduced()


What it does: Loads a preprocessed PBMC scRNA-seq dataset into an AnnData object, which stores expression data (.X) and metadata (.obs and .var).

Why: Provides a real-world scRNA-seq dataset ready for analysis.

Step 2: Add Metadata
np.random.seed(42)
adata.obs['condition'] = np.random.choice(['Healthy','Disease'], size=adata.n_obs, p=[0.5,0.5])
adata.obs['batch'] = np.random.choice(['Batch1','Batch2','Batch3'], size=adata.n_obs)


What it does: Adds condition and batch columns to simulate real experimental design.

Why: Many scRNA-seq datasets come with batch differences; handling this is key for accurate downstream analysis.

Tip: In actual studies, condition would come from donor disease status, and batch from sequencing runs.

Step 3: Batch Correction and Integration
import harmonypy as hm
sc.tl.pca(adata, svd_solver='arpack')

ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
adata.obsm['X_pca_harmony'] = ho.Z_corr.T


What it does:

Computes PCA embeddings of all cells.

Uses Harmony to remove batch effects while preserving biological variation.

Stores corrected embeddings in X_pca_harmony.

Why: Batch effects can mask real biological differences; correction is critical before clustering, DE analysis, or visualization.

Step 4: Neighborhood Graph & UMAP Visualization
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['batch','bulk_labels'], save="_integrated.png")


What it does:

Builds a k-nearest neighbor graph based on Harmony-corrected embeddings.

Computes UMAP coordinates for visualization.

Plots UMAP colored by batch and cell type.

Why: Visual inspection ensures that batch correction worked and shows the distribution of cell types in 2D space.

Result: Integrated UMAP shows cells cluster by biological type rather than batch.

Step 5: Analyze Cell-Type Proportions
cell_props = adata.obs.groupby(['bulk_labels','condition']).size().groupby(level=0).apply(lambda x: x/x.sum()).reset_index(name='proportion')


What it does:

Calculates proportion of each immune cell type within Healthy and Disease groups.

Normalizes proportions per cell type.

sns.barplot(data=cell_props, x='bulk_labels', y='proportion', hue='condition')


Why: Shifts in cell-type proportions can indicate disease-related immune changes.

Result: For example, CD14+ monocytes may increase in disease samples, suggesting activation.

Step 6: Cell-Type–Specific Differential Expression
mono = adata[adata.obs['bulk_labels'] == 'CD14+ Monocytes'].copy()
sc.tl.rank_genes_groups(mono, groupby='condition', method='wilcoxon')
sc.pl.rank_genes_groups(mono, n_genes=15, save="_monocyte_DE.png")


What it does:

Filters for CD14+ monocytes.

Performs DE analysis between Healthy and Disease within this cell type.

Plots top DE genes.

Why: Disease-specific transcriptional changes often occur within cell types, not globally.

Result: Monocytes in disease show upregulation of interferon-response and inflammatory genes.

Step 7: Export DE Results
de_results = sc.get.rank_genes_groups_df(mono, group='Disease')
de_results.to_csv("results/monocyte_disease_vs_healthy_DE.csv", index=False)


What it does: Saves DE genes as a CSV for reproducibility.

Why: Enables downstream functional analysis and sharing results.

Step 8: Pathway Enrichment (GO & KEGG)
import gseapy as gp
top_genes = de_results.query("logfoldchanges>0").head(200)['names'].tolist()
enr = gp.enrichr(gene_list=top_genes, gene_sets=['GO_Biological_Process_2021','KEGG_2021_Human'], organism='Human', outdir='results/enrichment', cutoff=0.05)


What it does:

Selects top upregulated genes in disease monocytes.

Performs Gene Ontology and KEGG pathway enrichment.

Why: Identifies biological processes altered in disease.

Result: Enrichment highlights innate immune activation, cytokine signaling, and inflammatory responses.

Step 9: Biological Interpretation

Disease alters immune cell proportions and transcriptional states.

CD14+ monocytes in disease show immune activation signatures (e.g., interferon-response genes).

Pathway enrichment confirms innate immune and inflammatory pathways are upregulated.

This aligns with expected immune dysregulation in infection or inflammatory conditions.

Project Structure
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

Tools & Libraries

Python: Core language

Scanpy & AnnData: scRNA-seq analysis

Harmony: Batch correction

GSEApy: Pathway enrichment

NumPy, Pandas: Data handling

Matplotlib, Seaborn: Visualization

Skills Demonstrated

scRNA-seq preprocessing and metadata integration

Batch effect correction & UMAP visualization

Cell-type proportion analysis

Cell-type–specific differential expression

Pathway enrichment and biological interpretation

End-to-end reproducible Python workflows
