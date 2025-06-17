# Differential gene expression analysis using limma-voom
# Author: Max Kluyev

# 1. Metadata formatting and outlier removal
filtered_metadata$stage <- factor(filtered_metadata$stage)
filtered_metadata$Braak <- factor(filtered_metadata$Braak, levels = 0:6)
filtered_metadata$sex <- factor(filtered_metadata$sex)
filtered_metadata$race <- factor(filtered_metadata$race)
filtered_metadata$ageDeath <- as.numeric(ifelse(as.character(filtered_metadata$ageDeath) == "90+", "91", as.character(filtered_metadata$ageDeath)))
filtered_metadata$pmi <- as.numeric(as.character(filtered_metadata$pmi))
filtered_metadata$RIN <- as.numeric(as.character(filtered_metadata$RIN))
filtered_metadata$sequencingBatch <- factor(filtered_metadata$sequencingBatch)
filtered_metadata$totalReads <- as.numeric(as.character(filtered_metadata$totalReads))

# Example: remove outlier samples
rownames_to_remove <- c("BM_10_572", "hB_RNA_13684", "hB_RNA_13444", "BM_10_731", "hB_RNA_12444")
filtered_metadata <- filtered_metadata[!rownames(filtered_metadata) %in% rownames_to_remove, ]
filtered_counts <- filtered_counts[, rownames(filtered_metadata)]

# 2. Design matrix and voom transformation
design <- model.matrix(~ 0 + Braak + sex + RIN + sequencingBatch, data = filtered_metadata)
dge <- DGEList(counts = filtered_counts)
dge <- calcNormFactors(dge)
voom_data <- voom(dge, design)
fit <- lmFit(voom_data, design)

# 3. Example: severe vs healthy (early_vs_0) contrast
contrasts <- makeContrasts(
  early_vs_0 = Braak6 - (Braak0 + Braak1 + Braak2 + Braak3 + Braak4)/5,
  levels = design
)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

results <- decideTests(fit2, method = "nestedF", adjust.method = "BH")
cat("Summary for severe vs healthy:\n")
print(summary(results))

top_genes <- topTable(fit2, coef = "early_vs_0", adjust.method = "BH", number = Inf)
significant_genes <- top_genes[top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 0.5, ]

# 4. All pairwise Braak contrasts
braak_levels <- levels(filtered_metadata$Braak)
contrast_list <- list()
for(i in 1:(length(braak_levels)-1)) {
  for(j in (i+1):length(braak_levels)) {
    contrast_name <- paste0("Braak", braak_levels[j], "_vs_Braak", braak_levels[i])
    contrast_expr <- paste0("Braak", braak_levels[j], " - Braak", braak_levels[i])
    contrast_list[[contrast_name]] <- contrast_expr
  }
}
contrasts <- do.call(makeContrasts, c(contrast_list, list(levels = design)))

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
results <- decideTests(fit2, method = "nestedF", adjust.method = "BH")

summary_table <- t(apply(results, 2, function(x) {
  c(Down = sum(x == -1), NotDE = sum(x == 0), Up = sum(x == 1))
}))
summary_table <- as.data.frame(summary_table)

top_genes_list <- list()
for(contrast_name in colnames(contrasts)) {
  top_genes_list[[contrast_name]] <- topTable(fit2, coef = contrast_name, adjust.method = "BH", number = Inf)
}
significant_genes_list <- lapply(top_genes_list, function(tab) {
  tab[tab$adj.P.Val < 0.05 & abs(tab$logFC) > 0.5, ]
})
significant_gene_counts <- sapply(significant_genes_list, nrow)
summary_table$Significant <- significant_gene_counts
cat("Summary for all Braak contrasts:\n")
print(summary_table)

# 5. Save significant genes for each contrast
for(contrast in names(significant_genes_list)) {
  df <- significant_genes_list[[contrast]]
  if(nrow(df) > 0) {
    df$Contrast <- contrast
    file_name <- paste0("Significant_genes_", contrast, ".csv")
    write.csv(df, file = file_name, row.names = TRUE)
  }
}

# 6. Null distribution for bootstrap FDR
residuals_values <- residuals(fit, y = voom_data$E)
null_voom <- voom_data
null_voom$E <- residuals_values
fit_null <- lmFit(null_voom, design)
contrasts <- makeContrasts(
  early_vs_0 = Braak6 - (Braak0 + Braak1 + Braak2 + Braak3 + Braak4)/5,
  levels = design
)
fit2_null <- contrasts.fit(fit_null, contrasts)
fit2_null <- eBayes(fit2_null)
test_result_null <- topTable(fit2_null, coef = "early_vs_0", adjust.method = "BH", number = Inf)
results_null <- decideTests(fit2_null, method="nestedF", adjust.method = "BH")
cat("Null model summary:\n")
print(summary(results_null))

# 7. P-value distributions
hist(top_genes$P.Value, main = "Real p-value distribution", xlab = "p-value")
hist(test_result_null$P.Value, main = "Null (residual) p-value distribution", xlab = "p-value")

# 8. Bootstrap FDR estimation
null_voom$targets <- cbind(
  null_voom$targets,
  filtered_metadata[colnames(null_voom$E), c("Braak", "sex", "RIN", "sequencingBatch")]
)
null_voom$design <- NULL

bootstrap_voom <- function(voom_obj, stage_col = "Braak") {
  if (!stage_col %in% colnames(voom_obj$targets)) {
    stop(paste("Column", stage_col, "not found in voom_obj$targets"))
  }
  stages <- voom_obj$targets[[stage_col]]
  categories <- unique(stages)
  boot_idx <- integer(0)
  for (cat in categories) {
    idx_cat <- which(stages == cat)
    n_cat <- length(idx_cat)
    boot_idx_cat <- sample(idx_cat, size = n_cat, replace = TRUE)
    boot_idx <- c(boot_idx, boot_idx_cat)
  }
  boot_voom <- voom_obj
  boot_voom$E <- voom_obj$E[, boot_idx, drop = FALSE]
  boot_voom$weights <- voom_obj$weights[, boot_idx, drop = FALSE]
  boot_voom$targets <- voom_obj$targets[boot_idx, , drop = FALSE]
  return(boot_voom)
}

# Bootstrap loop
p_values_list_boot <- list()
for (i in 1:500) {
  print(i)
  boot_voom <- bootstrap_voom(null_voom, stage_col = "Braak")
  boot_voom$targets$stage <- factor(filtered_metadata$Braak)
  boot_voom$targets$sex <- factor(boot_voom$targets$sex)
  boot_voom$targets$RIN <- as.numeric(as.character(boot_voom$targets$RIN))
  boot_voom$targets$sequencingBatch <- factor(boot_voom$targets$sequencingBatch)
  design <- model.matrix(~ 0 + Braak, data = boot_voom$targets, na.action = na.pass)
  colnames(design) <- make.names(colnames(design))
  fit <- lmFit(boot_voom, design)
  contrasts <- makeContrasts(
    early_vs_0 = Braak6 - (Braak0 + Braak1 + Braak2 + Braak3 + Braak4)/5,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  top_genes_boot <- topTable(fit2, number = Inf, adjust.method = "none", sort.by = "none")
  p_values_list_boot[[i]] <- top_genes_boot$P.Value
}

# Plot p-value histogram for bootstrap samples
library(ggplot2)
data <- data.frame(p_value = unlist(p_values_list_boot, use.names = FALSE))
plot <- ggplot(data, aes(x = p_value)) +
  geom_histogram(bins = 50, fill = 'blue', color = 'black', alpha = 0.7) +
  geom_vline(xintercept = 0.05, color = 'red', linetype = "dashed", linewidth = 1) +
  labs(title = "P-value distribution (bootstrap null)", x = "P-value", y = "Count") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
print(plot)

# 9. Bootstrap summary statistics
results_summary_bootstrap <- data.frame(
  Bootstrap = seq_along(p_values_list_boot),
  MinPValue = sapply(p_values_list_boot, function(p_vals) min(p_vals, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
print(results_summary_bootstrap)

# Count of real p-values greater than min in bootstrap
real_p_values <- top_genes$P.Value
results_summary_bootstrap$PValueGreaterThanMin <- sapply(
  results_summary_bootstrap$MinPValue,
  function(min_p_val) sum(real_p_values > min_p_val, na.rm = TRUE)
)

# Estimate pi0 (proportion of true nulls)
pi_zero <- median(results_summary_bootstrap$PValueGreaterThanMin, na.rm = TRUE) / nrow(top_genes)

# 10. Iterative determination of optimal p-value cutoff
q_target <- 0.05
tolerance <- 0.00001
max_iter <- 1000
p_lower <- min(top_genes$P.Value, na.rm = TRUE)
p_upper <- max(top_genes$P.Value, na.rm = TRUE)
iteration_count <- 0
best_p_val <- NA
best_q_diff <- Inf
best_q_value <- NA

while ((p_upper - p_lower) > tolerance && iteration_count < max_iter) {
  iteration_count <- iteration_count + 1
  p_mid <- (p_lower + p_upper) / 2
  count_significant_genes <- sum(top_genes$P.Value < p_mid, na.rm = TRUE)
  if (count_significant_genes == 0) {
    warning("No significant genes at this cutoff. Stopping search.")
    break
  }
  SignificantGenesBoot <- mean(sapply(p_values_list_boot, function(p_vals) {
    sum(p_vals < p_mid, na.rm = TRUE)
  }))
  q_value <- (SignificantGenesBoot / count_significant_genes) * pi_zero
  q_diff <- abs(q_value - q_target)
  cat("Iteration:", iteration_count,
      "p_mid:", p_mid,
      "q_value:", q_value,
      "n_significant:", count_significant_genes, "\n")
  if (q_diff < best_q_diff) {
    best_q_diff <- q_diff
    best_p_val <- p_mid
    best_q_value <- q_value
  }
  if (q_value > q_target) {
    p_upper <- p_mid
  } else {
    p_lower <- p_mid
  }
}
cat("\nBest solution found:\n")
cat("p-value cutoff:", best_p_val, "\n")
cat("q-value:", best_q_value, "\n")
cat("Number of significant genes:", sum(top_genes$P.Value < best_p_val, na.rm = TRUE), "\n")

# 11. Extract and filter significant genes by bootstrapped q-value
significant_results_boot <- top_genes[top_genes$P.Value < best_p_val, ]
filtered_genes_boot <- significant_results_boot[abs(significant_results_boot$logFC) > 0.5, ]
num_significant_genes_boot <- nrow(significant_results_boot)
