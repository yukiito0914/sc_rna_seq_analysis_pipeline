# scRNA-seq Analysis Pipeline
This repository contains an R-based single-cell RNA-seq (scRNA-seq) analysis pipeline implemented in R Markdown, with core data processing and clustering performed using the Seurat R package. The workflow covers data loading, quality control, filtering, normalization, feature selection, dimensionality reduction, clustering, integration, cell-type annotation, pseudotime inference, and downstream analyses such as differential expression and cell proportion comparison. The pipeline was developed and tested on hepatoblastoma datasets from Bondoc *et al.* (2021).

## Table of Contents

- [scRNA-seq Analysis Pipeline](#scrna-seq-analysis-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Data Organization](#data-organization)
  - [Running the Pipeline](#running-the-pipeline)
  - [Pipeline Overview](#pipeline-overview)
    - [1. Setup and Dependencies](#1-setup-and-dependencies)
    - [2. Setup and Dependencies](#2-setup-and-dependencies)
    - [3. Quality Control](#3-quality-control)
    - [4. Filtering and Doublet Detection](#4-filtering-and-doublet-detection)
    - [5. Normalization and Feature Selection](#5-normalization-and-feature-selection)
    - [6. Dimensionality Reduction](#6-dimensionality-reduction)
    - [7. Clustering and UMAP](#7-clustering-and-umap)
    - [8. Batch Integration (Harmony)](#8-batch-integration-harmony)
    - [9. Marker Gene and Annotation](#9-marker-gene-and-annotation)
    - [10. Pseudotime Analysis](#10-pseudotime-analysis)
    - [11. Cell Proportion Analysis](#11-cell-proportion-analysis)
    - [12. Differential Expression Example](#12-differential-expression-example)
  - [Outputs](#outputs)
  - [Reference](#reference)

## Prerequisites

- **R** (>= 4.0)  
- **RStudio** (recommended for R Markdown execution)  
- Linux or macOS environment with multicore support  

## Installation

Install the required R packages:

```r
install.packages(c(
  "Seurat", "scDblFinder", "ggplot2", "harmony", "dplyr", "future", "patchwork", "pheatmap"
))
BiocManager::install(c(
  "SingleR", "scran", "slingshot", "CellChat", "celldex"
))
```

## Data Organization

Place the 10X Genomics filtered feature-barcode matrices under:

```r
refs/
├── HB17_background_filtered_feature_bc_matrix/
├── HB17_PDX_filtered_feature_bc_matrix/
├── HB17_tumor_filtered_feature_bc_matrix/
├── HB30_PDX_filtered_feature_bc_matrix/
├── HB30_tumor_filtered_feature_bc_matrix/
├── HB53_background_filtered_feature_bc_matrix/
└── HB53_tumor_filtered_feature_bc_matrix/
```

## Running the Pipeline

Open and knit the ```seurat.Rmd``` file in RStudio, or run:

```r
Rscript -e "rmarkdown::render('final_project_sc_rna_analysis.Rmd')"
```

This will produce an HTML report with embedded code, figures, and results.

## Pipeline Overview

### 1. Setup and Dependencies
- Load libraries: Seurat, scDblFinder, ggplot2, harmony, dplyr, SingleR, future, slingshot, CellChat, patchwork, pheatmap.
- Configure multicore support with ```plan("multicore", workers = 8)```.

### 2. Setup and Dependencies
- Read 10X data for background liver, tumor, and PDX samples using ```Read10X``` and ```CreateSeuratObject```.
- Assign sample metadata.

### 3. Quality Control
- Compute mitochondrial percentage (```percent.mt```) per cell.
- Visualize QC metrics (```nFeature_RNA```, ```nCount_RNA```, ```percent.mt```) via violin plots and scatter plots.

### 4. Filtering and Doublet Detection
- Subset cells based on feature, count, and mitochondrial thresholds.
- Detect and remove doublets with ```scDblFinder```.
- Summarize cell and gene counts before/after filtering.

### 5. Normalization and Feature Selection
- Merge filtered objects into one Seurat object.
- Normalize data with ```NormalizeData``` (LogNormalize) or optionally ```SCTransform```.
- Identify 2,000 highly variable genes (```FindVariableFeatures```).

### 6. Dimensionality Reduction
- Scale all genes with ```ScaleData```.
- Perform PCA (```RunPCA```) and select appropriate PCs using ```ElbowPlot```.

### 7. Clustering and UMAP
- Build a kNN graph (```FindNeighbors```) on top PCs.
- Identify clusters (```FindClusters```, resolution = 0.2).
- Compute UMAP embedding (```RunUMAP```) and visualize clusters and sample origin.

### 8. Batch Integration (Harmony)
- Correct batch effects across samples using ```RunHarmony```.
- Re-run neighbor finding, clustering, and UMAP on Harmony embeddings.

### 9. Marker Gene and Annotation
- Find cluster markers (```FindAllMarkers```).
- Automatic annotation with ```SingleR``` against HumanPrimaryCellAtlas.
- Manual cluster labeling via canonical marker heatmaps and violin plots.

### 10. Pseudotime Analysis
- Convert to SingleCellExperiment.
- Infer trajectories with ```slingshot``` on UMAP.
- Visualize pseudotime gradients.

### 11. Cell Proportion Analysis
- Compare cell-type composition between tumor and PDX samples.
- Barplots of cluster percentages by sample group.

### 12. Differential Expression Example
- Subset CD8⁺ T cells and perform DE analysis (```FindMarkers```) between tumor and PDX.
- Generate volcano plot and heatmap of top DE genes.

## Outputs
- HTML report: ```final_project_sc_rna_analysis.html```
- Marker heatmap PNG: ```heatmap_top5_markers.png```
- Figures for QC, UMAP, integration, pseudotime, DE analyses.

## Reference
Bondoc, A., Glaser, K., Jin, K., Lake, C., Cairo, S., Geller, J., ... & Aronow, B. (2021). Identification of distinct tumor cell populations and key genetic mechanisms through single cell sequencing in hepatoblastoma. *Communications Biology, 4*(1), 1049.
