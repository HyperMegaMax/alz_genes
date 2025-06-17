# RNA-seq Data Analysis for Alzheimer's Disease Research

## Author
Maxim Kluyev

## Project Description
This repository contains R scripts and analysis results for the RNA-seq data from the MSBB dataset, aimed at investigating long-term gene expression changes associated with Alzheimer's Disease progression.

## Contents
- Data preprocessing and normalization
- Differential gene expression analysis (limma-voom)
- Bootstrap-based false discovery rate (FDR) correction
- PCA and variance partitioning analysis
- Visualization of analytical results

## Repository Structure
- **data/** – raw and processed data.
- **scripts/** – R scripts for analysis.
- **results/** – output tables from analyses.
- **figures/** – plots and visualizations.
- **docs/** – supplementary documentation and reports.

## Reproduction Instructions
1. Clone the repository.
2. Install the required packages using the provided script (`scripts/setup.R`).
3. Run the analysis scripts in the specified order from the `scripts` directory.

## Dependencies
Analysis is conducted using R (version 4.x) with libraries from Bioconductor and CRAN:
- limma, edgeR, DESeq2, variancePartition
- ggplot2, pheatmap, FactoMineR, factoextra

## License
The project is available under the MIT License—see [LICENSE](LICENSE).

## Citation
If you use this code or its results, please cite as follows:

> Maxim Kluyev. (2025). *Differential Gene Expression in Alzheimer's Disease* [GitHub Repository]. URL: <repository-url>

