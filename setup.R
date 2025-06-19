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
