# RNA-seq Data Analysis for Alzheimer's Disease Research

## Author
Maxim Kluyev

## Project Description
This repository contains R scripts and analysis results for the RNA-seq data from the MSBB dataset, aimed at investigating long-term gene expression changes associated with Alzheimer's Disease progression.

## Contents
- Data preprocessing and normalization
- PCA and variance partitioning analysis
- Differential gene expression analysis (limma-voom)
- Bootstrap-based false discovery rate (FDR) correction
- Network analysis
- FDR correction methods comparisson

## Data Requirements

This project uses post‑mortem RNA‑seq data from the Mount Sinai Brain Bank (MSBB) and the Religious Orders Study & Memory and Aging Project (ROSMAP).  
**These datasets must be obtained separately before analysis.** Both are publicly available through the AMP‑AD Knowledge Portal (Synapse):

- **MSBB**: available at [Synapse syn3159438](https://www.synapse.org/Synapse:syn3159438)
- **ROSMAP**: available at [Synapse syn3219045](https://www.synapse.org/Synapse:syn3219045)

Refer to the listed Synapse IDs and the AMP‑AD documentation to request and download the necessary clinical, metadata, and expression files.


## Repository Structure
- **scripts/** – Directory containing main R scripts for analysis:
    - `preprocessing.R` – data import, cleaning, and sample filtering
    - `eda_pca_outliers.R` – exploratory data analysis, PCA, and outlier detection
    - `differential_expression_limma_boot.R` – differential expression analysis using limma with bootstrap-based FDR correction
    - `FDR_correction_testing.R` – multiple testing correction and benchmarking (including BH, BY, Storey, bootstrap)
    - `network_analysis.R` – construction and clustering of gene co-expression or interaction networks
- **results/** – output tables from analyses.
    - `new_clusters_xl2.xlsx` – Results of network analysis and Louvain clustering.
    - `very_important_genes_1004_3.xlsx` – Results of differential expression analysis.

## Reproduction Instructions
1. Clone the repository.
2. Install the required packages using the provided script (`scripts/setup.R`).
3. Run the analysis scripts in the specified order from the `scripts` directory.


## Dependencies

This analysis is conducted using **R** (version 4.3.3) with packages from both **Bioconductor** and **CRAN**.

All code assumes the following packages are installed and up-to-date:

### Core Differential Expression & Preprocessing
- `limma` – linear modeling and voom transformation
- `edgeR` – count normalization (TMM) and DGE analysis
- `DESeq2` – alternative DGE pipeline
- `variancePartition` – variance decomposition and covariate assessment
- `NOISeq` – (optional) for additional normalization/comparison

### Data Wrangling & Utilities
- `dplyr`, `tidyr`, `tibble`, `readr` – fast data import, manipulation, and reshaping
- `reshape2` – data reshaping (legacy)
- `stringr` – string utilities

### Visualization
- `ggplot2` – publication-quality plots
- `pheatmap` – heatmaps
- `FactoMineR`, `factoextra` – PCA and exploratory analysis

### Multiple Testing & FDR Correction
- `qvalue` – Storey's FDR/q-value estimation

### Synthetic Data & Simulation
- `compcodeR` – simulation of RNA-seq count data for benchmarking

### Network Analysis & Graphs
- `igraph` – network creation, metrics, clustering (Louvain, degree)
- `ggraph`, `tidygraph` – advanced network visualization

### Cell Type Annotation
- `clustermole` (Bioconductor) – annotation of gene clusters by cell type

### Package Management
- `BiocManager` – for installing Bioconductor packages

---

### Installation

To install all required packages, run in R:

```r
# Install BiocManager if not already present
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor packages
BiocManager::install(c(
  "limma", "edgeR", "DESeq2", "variancePartition", 
  "NOISeq", "clustermole", "compcodeR"
))

# CRAN packages
install.packages(c(
  "dplyr", "tidyr", "tibble", "readr", "reshape2", "stringr",
  "ggplot2", "pheatmap", "FactoMineR", "factoextra",
  "igraph", "ggraph", "tidygraph", "qvalue"
))

```

## License
The project is available under the MIT License—see [LICENSE](LICENSE).

## Citation
If you use this code or its results, please cite as follows:

> Maxim Kluyev. (2025). *Differential Gene Expression in Alzheimer's Disease* [GitHub Repository]. URL: <repository-url>

