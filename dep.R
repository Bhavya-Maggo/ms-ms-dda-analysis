# 1. Load libraries----------------------------------------------------
library(DEP)
library(tidyverse)
library(SummarizedExperiment)
library(vsn)
library (ggplot2)
library(pheatmap)

# 2. Read data------------------------------------------------------------
data <- read.delim("proteinGroups.txt", check.names = FALSE)

data <- filter(data,
                      Reverse != "+",
                      `Potential contaminant` != "+")

experimental_design <- read.delim("experiment_annotation.txt", check.names = FALSE)


# 3. Handle duplicate gene names------------------------------------------------
data_unique <- make_unique(data, "Gene names", "Protein IDs", delim = ";")

# 4. Make SummarizedExperiment--------------------------------------------
LFQ_columns <- grep("LFQ.", colnames(data_unique))

data_se <- make_se(data_unique, LFQ_columns, experimental_design)

# 5. QC plots-------------------------------------------------------------
plot_frequency(data_se)   # how many samples each protein appears in
plot_numbers(data_se)     # proteins per sample

# 6. Filter missing values--------------------------------------------------
data_filt <- filter_missval(data_se, thr = 1)
meanSdPlot(assay(data_filt))
plot_numbers(data_filt)
plot_coverage(data_filt)

# 7. Normalize------------------------------------------------------------
data_norm <- normalize_vsn(data_filt)
meanSdPlot(assay(data_norm))
plot_normalization(data_filt, data_norm)

# 8. Impute-------------------------------------------------------------
plot_missval(data_filt) 
plot_detect(data_filt)  

data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

# 9. Statistical testing--------------------------------------------------- 
# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Control")

# Add significance flags (FDR 0.05, fold change 1.5)
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

table(rowData(dep)$significant)# How many significant per contrast?

# 10. Results table-------------------------------------------------
res <- get_results(dep)   
colnames(res)           

res %>% filter(significant) %>% nrow()

# 11. Visualization of the results----------------------------------------
## PCA plot
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

## Correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

## Heatmap: Centered: shows expression relative to mean per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

## Heatmap: Contrast: shows fold changes between conditions directly
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

## Volcano plots
contrasts_to_plot <- c(
  "X.HS276_vs_Control",
  "X.LLOMe5m_vs_Control",
  "LLOMe30m_vs_Control",
  "X.HS276_LLOMe5m_vs_Control",
  "X.HS276_LLOMe30m_vs_Control"
)
for (ct in contrasts_to_plot) {
  
  ratio_col <- paste0(ct, "_ratio")    # log2 fold change column
  pval_col  <- paste0(ct, "_p.val")   # p-value column
  sig_col   <- paste0(ct, "_significant")
  
  # Check the column exists before plotting
  if (!ratio_col %in% colnames(res)) {
    message("Column not found: ", ratio_col, " — check colnames(res)")
    next
  }
  
  volcano_data <- res %>%
    mutate(
      log2FC = .data[[ratio_col]],
      pval   = .data[[pval_col]],
      sig    = .data[[sig_col]]
    ) %>%
    filter(!is.na(log2FC), !is.na(pval))
  
  p <- ggplot(volcano_data, aes(x = log2FC, y = -log10(pval))) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
               linetype = "dashed", color = "#EF9F27", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "#EF9F27", linewidth = 0.5) +
    scale_color_manual(values = c("FALSE" = "#B4B2A9", "TRUE" = "#D85A30"),
                       labels = c("Not significant", "Significant"),
                       name = NULL) +
    # Label top 10 significant hits
    ggrepel::geom_text_repel(
      data = . %>% filter(sig == TRUE) %>% slice_min(pval, n = 10),
      aes(label = name), size = 2.8, max.overlaps = 15
    ) +
    labs(title = gsub("_", " ", ct),
         x = "log2 fold change",
         y = "-log10 p-value") +
    theme_bw(base_size = 12)
  
  print(p)
  ggsave(paste0("volcano_", ct, ".pdf"), p, width = 7, height = 6)
}

# 12. Condition overlap
plot_cond(dep)
