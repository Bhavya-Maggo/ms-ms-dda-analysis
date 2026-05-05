##########################################################################
# Load libraries
library(limma)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(dplyr)
library (ggrepel)

##########################################################################
# Read the file
mq <- read.delim("proteinGroups.txt", sep = "\t", header = TRUE)

# Remove contaminants, reverse hits, and only-by-site proteins
mq <- mq[mq$Potential.contaminant != "+", ]
mq <- mq[mq$Reverse != "+", ]
mq <- mq[mq$Only.identified.by.site != "+", ]

clean_names <- function(x) {
  x <- sapply(strsplit(as.character(x), ";"), `[`, 1)
  idx <- is.na(x) | x == ""
  x[idx] <- paste0("Protein_", seq_len(sum(idx)))
  
  x
}

# Extract LFQ intensity columns
lfq_cols <- grep("LFQ.intensity.", colnames(mq), value = TRUE)
mat <- as.matrix(mq[, lfq_cols])
rownames(mat) <- make.unique(clean_names(mq$Gene.names))  # make.unique catches any remaining duplicates

##########################################################################
# FILTER → NORMALIZE → IMPUTE
mat[mat == 0] <- NA
mat_log <- log2(mat)

group <- factor(c(
  "Control",        "Control",        "Control",
  "HS276",          "HS276",          "HS276",
  "HS276_LLOMe30m", "HS276_LLOMe30m", "HS276_LLOMe30m",
  "HS276_LLOMe5m",  "HS276_LLOMe5m",  "HS276_LLOMe5m",
  "LLOMe30m",       "LLOMe30m",       "LLOMe30m",
  "LLOMe5m",        "LLOMe5m",        "LLOMe5m"
))

# Filter: ≥ 2 valid values in at least one group
keep <- sapply(levels(group), function(g) {
  cols <- which(group == g)
  rowSums(!is.na(mat_log[, cols])) >= 2
})
mat_filt <- mat_log[rowSums(keep) >= 1, ]

# Normalize
mat_norm <- normalizeBetweenArrays(mat_filt, method = "quantile")
boxplot(mat_norm, las = 2, main = "Post-normalization", cex.axis = 0.6)

# Impute
set.seed(42)
mat_imp <- mat_norm
for (col in 1:ncol(mat_imp)) {
  nas <- is.na(mat_imp[, col])
  if (sum(nas) == 0) next
  col_min <- min(mat_imp[, col], na.rm = TRUE)
  col_sd  <- sd(mat_imp[, col],  na.rm = TRUE)
  mat_imp[nas, col] <- rnorm(sum(nas),
                             mean = col_min - 1.8 * col_sd,
                             sd   = 0.3   * col_sd)
}

##########################################################################
# LIMMA MODEL
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(mat_imp, design)   # <-- mat_imp not mat_norm

contrast_mat <- makeContrasts(
  HS276_vs_Control          = HS276 - Control,
  LLOMe5m_vs_Control        = LLOMe5m - Control,
  LLOMe30m_vs_Control       = LLOMe30m - Control,
  HS276_LLOMe5m_vs_Control  = HS276_LLOMe5m - Control,
  HS276_LLOMe30m_vs_Control = HS276_LLOMe30m - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2, trend = TRUE)

# Save CSVs
for (contrast in colnames(contrast_mat)) {
  res <- topTable(fit2, coef = contrast, number = Inf, adjust.method = "BH")
  write.csv(res, paste0("results_", contrast, ".csv"))
}

##########################################################################
# Volcano plots
for (contrast in colnames(contrast_mat)) {
  res <- topTable(fit2, coef = contrast, number = Inf, adjust.method = "BH")
  res$gene <- rownames(res)
  res$adj.P.Val[res$adj.P.Val == 0] <- 1e-300
  
  res$Regulation <- "Not significant"
  res$Regulation[res$adj.P.Val < 0.05 & res$logFC >  1] <- "Upregulated"
  res$Regulation[res$adj.P.Val < 0.05 & res$logFC < -1] <- "Downregulated"
  
  top_labels <- res %>%
    filter(Regulation != "Not significant") %>%
    slice_min(adj.P.Val, n = 15)
  
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
    geom_point(size = 1.8, alpha = 0.7) +
    geom_text_repel(data = top_labels, aes(label = gene),
                    size = 2.8, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(
      "Upregulated"     = "#D85A30",
      "Downregulated"   = "#378ADD",
      "Not significant" = "#B4B2A9")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    labs(title = gsub("_", " ", contrast),
         x = "log2 Fold Change", y = "-log10 adjusted p-value") +
    theme_bw(base_size = 13)
  
  ggsave(paste0("Volcano_", contrast, ".png"), p, width = 7, height = 6, dpi = 300)
}

##########################################################################
# Heatmaps
anno_col <- data.frame(Group = group)
rownames(anno_col) <- colnames(mat_imp)

anno_colors <- list(Group = c(
  "Control"        = "#888780",
  "HS276"          = "#7F77DD",
  "HS276_LLOMe30m" = "#1D9E75",
  "HS276_LLOMe5m"  = "#5DCAA5",
  "LLOMe30m"       = "#D85A30",
  "LLOMe5m"        = "#F0997B"
))

for (contrast in colnames(contrast_mat)) {
  res <- topTable(fit2, coef = contrast, number = Inf, adjust.method = "BH")
  res$gene <- rownames(res)
  
  top_genes <- res %>%
    filter(adj.P.Val < 0.05) %>%       # only significant
    slice_min(adj.P.Val, n = 50) %>%
    pull(gene)
  
  mat_sig <- mat_imp[rownames(mat_imp) %in% top_genes, ]
  
  if (nrow(mat_sig) < 2) {
    message("Skipping ", contrast, " — fewer than 2 significant proteins")
    next
  }
  
  mat_scaled <- t(scale(t(mat_sig)))
  
  pheatmap(mat_scaled,
           annotation_col  = anno_col,
           annotation_colors = anno_colors,
           show_rownames   = nrow(mat_scaled) <= 40,  # show names if not too crowded
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           color = colorRampPalette(c("#378ADD", "white", "#D85A30"))(100),
           border_color    = NA,
           main            = gsub("_", " ", contrast),
           filename        = paste0("Heatmap_", contrast, ".png"),
           width = 7, height = 8)
}

##########################################################################
# PCA
ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = Sample), size = 2.8, show.legend = FALSE) +
  
  stat_ellipse(
    data = dplyr::filter(pca_df, Group %in% valid_groups),
    aes(x = PC1, y = PC2, group = Group),
    linetype = 2,
    linewidth = 0.4
  ) +
  
  scale_color_manual(values = c(
    "Control"        = "#888780",
    "HS276"          = "#7F77DD",
    "HS276_LLOMe30m" = "#1D9E75",
    "HS276_LLOMe5m"  = "#5DCAA5",
    "LLOMe30m"       = "#D85A30",
    "LLOMe5m"        = "#F0997B"
  )) +
  labs(
    title = "PCA — all 18 samples",
    x = paste0("PC1 (", var_ex[1], "%)"),
    y = paste0("PC2 (", var_ex[2], "%)")
  ) +
  theme_bw(base_size = 13)
