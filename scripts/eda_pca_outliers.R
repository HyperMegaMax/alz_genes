# eda_pca_outliers.R
# Exploratory data analysis, PCA, and outlier detection for MSBB RNA-seq data

# ---------------------- #
# 1. Normalization (TMM + voom)
# ---------------------- #

# Create DGEList object and TMM normalization
dge <- DGEList(counts = filtered_counts)
dge <- calcNormFactors(dge)

# voom transformation (for downstream analysis)
voom_data <- voom(dge, plot = TRUE)
voom_counts <- voom_data$E
filtered_counts <- voom_counts

# ---------------------- #
# 2. Metadata formatting
# ---------------------- #
filtered_metadata$stage <- factor(filtered_metadata$stage, levels = c("healthy", "mild", "moderate", "severe"))
filtered_metadata$sex   <- factor(filtered_metadata$sex)
filtered_metadata$race  <- factor(filtered_metadata$race)
filtered_metadata$ageDeath <- as.numeric(ifelse(as.character(filtered_metadata$ageDeath) == "90+", "91", as.character(filtered_metadata$ageDeath)))
filtered_metadata$pmi   <- as.numeric(as.character(filtered_metadata$pmi))
filtered_metadata$RIN   <- as.numeric(as.character(filtered_metadata$RIN))
filtered_metadata$sequencingBatch <- factor(filtered_metadata$sequencingBatch)
filtered_metadata$totalReads      <- as.numeric(as.character(filtered_metadata$totalReads))

# ---------------------- #
# 3. Covariate variance analysis (variancePartition)
# ---------------------- #
form <- ~ stage + sex + race + ageDeath + pmi + RIN + sequencingBatch + totalReads
vp <- fitExtractVarPartModel(voom_data$E, form, filtered_metadata)
plotVarPart(vp)

# ---------------------- #
# 4. Prepare numeric covariates for PCA correlations
# ---------------------- #
library(dplyr)
metadata_numeric <- filtered_metadata %>%
  mutate(
    ageDeath = as.numeric(as.character(ageDeath)),
    pmi      = as.numeric(as.character(pmi)),
    RIN      = as.numeric(as.character(RIN)),
    totalReads = as.numeric(as.character(totalReads)),
    stage = as.numeric(as.factor(stage)),
    sex   = as.numeric(as.factor(sex)),
    race  = as.numeric(as.factor(race)),
    sequencingBatch = as.numeric(as.factor(sequencingBatch))
  )

# ---------------------- #
# 5. PCA and correlation with covariates
# ---------------------- #
library(pheatmap)
expr <- filtered_counts

# PCA: rows = genes, columns = samples (transpose matrix)
pca_res <- prcomp(t(expr), scale. = TRUE)
pc_scores <- as.data.frame(pca_res$x)[, 1:10]
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
pc_names <- paste0("PC", 1:10, " (", round(var_explained[1:10], 1), "%)")

cov_names <- c("stage", "sex", "race", "ageDeath", "pmi", "RIN", "sequencingBatch", "totalReads")
n_pcs <- 10; n_cov <- length(cov_names)

cor_matrix <- matrix(NA, nrow = n_pcs, ncol = n_cov, dimnames = list(pc_names, cov_names))
pval_matrix <- matrix(NA, nrow = n_pcs, ncol = n_cov, dimnames = list(pc_names, cov_names))
common_samples <- intersect(rownames(pc_scores), rownames(metadata_numeric))

for (i in 1:n_pcs) {
  for (j in 1:n_cov) {
    cov_vec <- metadata_numeric[common_samples, cov_names[j]]
    if (!is.numeric(cov_vec)) cov_vec <- suppressWarnings(as.numeric(as.character(cov_vec)))
    pc_vec  <- pc_scores[common_samples, i]
    complete_idx <- complete.cases(pc_vec, cov_vec)
    if (sum(complete_idx) >= 3) {
      res <- cor.test(pc_vec[complete_idx], cov_vec[complete_idx], method = "pearson")
      cor_matrix[i, j] <- res$estimate
      pval_matrix[i, j] <- res$p.value
    }
  }
}
# FDR correction
p_adjusted <- matrix(p.adjust(as.vector(pval_matrix), method = "BH"),
                     nrow = n_pcs, ncol = n_cov, dimnames = list(pc_names, cov_names))
sig_matrix <- ifelse(p_adjusted < 0.1, "*", "")
display_matrix <- matrix(paste0(round(cor_matrix, 2), sig_matrix),
                         nrow = n_pcs, ncol = n_cov, dimnames = list(pc_names, cov_names))
finite_vals <- cor_matrix[is.finite(cor_matrix)]
my_breaks <- seq(min(finite_vals), max(finite_vals), length.out = 101)
pheatmap(cor_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Pearson Correlation: 10 PCs vs Covariates",
         display_numbers = display_matrix,
         number_format = "%.2f",
         fontsize_number = 8,
         breaks = my_breaks)

# ---------------------- #
# 6. Detailed PCA for outlier search
# ---------------------- #
library(ggplot2)
library(tibble)
library(dplyr)

# Save metadata for visualization
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$SampleID <- rownames(pca_df)
pca_df <- left_join(pca_df, filtered_metadata %>% rownames_to_column("SampleID"), by = "SampleID")

# Euclidean distance to center
center_PC1 <- mean(pca_df$PC1)
center_PC2 <- mean(pca_df$PC2)
pca_df <- pca_df %>%
  mutate(dist_to_center = sqrt((PC1 - center_PC1)^2 + (PC2 - center_PC2)^2)) %>%
  arrange(desc(dist_to_center))
n_outliers <- 10
outlier_samples <- pca_df$SampleID[1:n_outliers]
pca_df <- pca_df %>%
  mutate(label = ifelse(SampleID %in% outlier_samples, SampleID, ""))

# PCA plot with outlier labels
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2,
                         color = sequencingBatch,
                         shape = stage,
                         size = RIN)) +
  geom_point() +
  geom_text(aes(label = label),
            hjust = 1.1, vjust = 1.1,
            size = 3, color = "black") +
  theme_bw() +
  labs(title = "PCA of normalized RNA-seq data with outliers labeled",
       x = "PC1", y = "PC2")
print(p1)

# ---------------------- #
# 7. Outlier detection by Mahalanobis distance
# ---------------------- #
data_matrix <- as.matrix(pca_df[, c("PC1", "PC2")])
center <- colMeans(data_matrix)
cov_matrix <- cov(data_matrix)
mahal_dist <- mahalanobis(x = data_matrix, center = center, cov = cov_matrix)
threshold <- qchisq(0.995, df = 2)
outlier_indices <- which(mahal_dist > threshold)
sorted_indices <- outlier_indices[order(mahal_dist[outlier_indices], decreasing = TRUE)]
outlier_samples_mahal <- pca_df$SampleID[sorted_indices]
cat("Outliers by Mahalanobis distance:", outlier_samples_mahal, "\n")

# Example: how to remove outliers
# rownames_to_remove <- outlier_samples_mahal[1:3]
# filtered_metadata <- filtered_metadata[!rownames(filtered_metadata) %in% rownames_to_remove, ]
# filtered_counts <- filtered_counts[, rownames(filtered_metadata)]

# End of EDA/PCA script
