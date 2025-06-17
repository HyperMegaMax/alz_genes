# preprocess_msbb.R
# Preprocessing script for MSBB RNA-seq data
# Author: Max Kluyev

# ------------------ #
# 1. Install and load required packages ----
# ------------------ #
packages <- c(
  "BiocManager", "Biobase", "limma", "ggplot2", "pheatmap", "reshape2",
  "FactoMineR", "factoextra", "dplyr", "tidyr", "tibble",
  "DESeq2", "edgeR", "variancePartition", "car"
)
bioc_packages <- c("Biobase", "limma", "DESeq2", "edgeR", "variancePartition", "polyester")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% bioc_packages) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Load libraries
library(Biobase)
library(limma)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(tidyr)
library(tibble)
library(variancePartition)
library(car)

# ------------------ #
# 2. Data Loading ----
# ------------------ #

# Set correct paths to your data folder!
meta_dir <- "data/raw/metadata/"
counts_file <- "data/raw/MSSM_FP_STG_PHG_IFG_Counts.tsv"

MSBB_individual_metadata   <- read.csv(paste0(meta_dir, "MSBB_individual_metadata.csv"))
MSBB_biospecimen_metadata <- read.csv(paste0(meta_dir, "MSBB_biospecimen_metadata.csv"))
MSBB_assay_rnaSeq_metadata <- read.csv(paste0(meta_dir, "MSBB_assay_rnaSeq_metadata_corrected.csv"))
count_data <- read.table(counts_file, sep = "\t", header = TRUE, check.names = FALSE)

# ------------------ #
# 3. Metadata Merging and Filtering ----
# ------------------ #

metadata <- MSBB_biospecimen_metadata %>%
  inner_join(MSBB_individual_metadata, by = "individualID") %>%
  inner_join(MSBB_assay_rnaSeq_metadata, by = "specimenID") %>%
  filter(assay == "rnaSeq")

# Remove samples without Braak score
metadata <- metadata[!is.na(metadata$Braak), ]

# Assign stage factor
classify_stage <- function(Braak_stage) {
  if (is.na(Braak_stage)) {
    return(NA)
  } else if (Braak_stage %in% 0:4) {
    return("healthy")
  } else if (Braak_stage == 5) {
    return("moderate")
  } else if (Braak_stage == 6) {
    return("severe")
  } else {
    return(NA)
  }
}
metadata$stage <- mapply(classify_stage, metadata$Braak)

rownames(metadata) <- metadata$specimenID

# Reformat counts
rownames(count_data) <- count_data$ensembl_gene_id
count_data$ensembl_gene_id <- NULL

# ------------------ #
# 4. Brodmann Area Selection (optional) ----
# ------------------ #
# Uncomment the area you want to analyze. By default, use Brodmann area 10.
metadata <- metadata[metadata$BrodmannArea == "10", ]

# ------------------ #
# 5. Synchronize samples and filter NAs ----
# ------------------ #
common_samples <- intersect(rownames(metadata), colnames(count_data))
filtered_metadata <- metadata[common_samples, , drop = FALSE]
filtered_counts <- count_data[, common_samples]

# Remove samples with missing covariates
covariates <- c("Braak", "sex", "race", "ageDeath", "pmi", "RIN", "sequencingBatch", "totalReads")
filtered_metadata <- filtered_metadata[complete.cases(filtered_metadata[, covariates]), ]
filtered_counts   <- filtered_counts[, rownames(filtered_metadata)]

filtered_metadata <- filtered_metadata[, covariates]  # Leave only relevant columns

# ------------------ #
# 6. Normalization (TMM + voom) ----
# ------------------ #
dge <- DGEList(counts = filtered_counts)
dge <- calcNormFactors(dge)
voom_data <- voom(dge, plot = TRUE)
voom_counts <- voom_data$E

# Save normalized data for downstream analysis
saveRDS(filtered_metadata, file = "data/processed/filtered_metadata.rds")
saveRDS(voom_counts, file = "data/processed/voom_counts.rds")

# ------------------ #
# 7. Covariate Formatting ----
# ------------------ #
filtered_metadata$stage <- factor(filtered_metadata$stage, levels = c("healthy", "mild", "moderate", "severe"))
filtered_metadata$sex <- factor(filtered_metadata$sex)
filtered_metadata$race <- factor(filtered_metadata$race)
filtered_metadata$ageDeath <- as.numeric(ifelse(as.character(filtered_metadata$ageDeath) == "90+", "91", as.character(filtered_metadata$ageDeath)))
filtered_metadata$pmi <- as.numeric(as.character(filtered_metadata$pmi))
filtered_metadata$RIN <- as.numeric(as.character(filtered_metadata$RIN))
filtered_metadata$sequencingBatch <- factor(filtered_metadata$sequencingBatch)
filtered_metadata$totalReads <- as.numeric(as.character(filtered_metadata$totalReads))

# Optionally, save cleaned and formatted metadata
saveRDS(filtered_metadata, file = "data/processed/filtered_metadata_formatted.rds")

# ------------------ #
# End of preprocessing script

