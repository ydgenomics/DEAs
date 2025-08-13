# Title: volcano_visual.R
# Date: 2025-05-14
# csv: avg_log2FC, p_val_adj, gene, cluster

library(optparse)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define command-line options
option_list <- list(
  make_option(c("-g", "--gene_csv"), type = "character", default = "/data/work/peanut/DEA/markers_peanut.csv", help = "Path to the gene CSV file"),
  make_option(c("-c", "--coef_col"), type = "character", default = "avg_log2FC", help = "Column name for coefficients"),
  make_option(c("-p", "--pval_col"), type = "character", default = "p_val_adj", help = "Column name for p-values"),
  make_option(c("-t", "--coef_threshold"), type = "numeric", default = 1, help = "Threshold for coefficients"),
  make_option(c("-q", "--pval_threshold"), type = "numeric", default = 0.05, help = "Threshold for p-values")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Assign parsed options to variables
gene_csv <- opt$gene_csv
coef_col <- opt$coef_col
pval_col <- opt$pval_col
coef_threshold <- opt$coef_threshold
pval_threshold <- opt$pval_threshold


markers <- read.csv(gene_csv, header = TRUE, stringsAsFactors = FALSE)
head(markers)
# csv: avg_log2FC, p_val_adj, gene, cluster

plot_volcano <- function(data, coef_col, pval_col, coef_threshold, pval_threshold, title, filename) {
  #Significance
  data <- data %>%
    mutate(Significance = case_when(
      (!!sym(coef_col) > coef_threshold & !!sym(pval_col) < pval_threshold) ~ "Significant Up",
      (!!sym(coef_col) < -coef_threshold & !!sym(pval_col) < pval_threshold) ~ "Significant Down",
      TRUE ~ "Not Significant"
    ))
  
  # deal zero p-value
  zero_pval_genes <- filter(data, !!sym(pval_col) == 0)
  if (nrow(zero_pval_genes) > 0) {
    min_nonzero_pval <- min(filter(data, !!sym(pval_col) > 0)[[pval_col]], na.rm = TRUE)
    if (is.na(min_nonzero_pval)) {
      min_nonzero_pval <- 1e-10
    }
    data <- mutate(data, !!sym(pval_col) := ifelse(!!sym(pval_col) == 0, min_nonzero_pval, !!sym(pval_col)))
  }
  
  # -log10(p-value)
  data <- mutate(data, neg_log_pval = -log(!!sym(pval_col), base = 10))
  
  # plot volcano
  ggplot(data, aes(x = !!sym(coef_col), y = neg_log_pval, color = Significance)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-coef_threshold, coef_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log(pval_threshold, base = 10), linetype = "dashed", color = "black") +
    labs(title = title, x = coef_col, y = "-Log10(P-value)", color = "Significance") +
    theme_minimal() +
    theme(legend.position = "right") +
    ggrepel::geom_text_repel(
      data = data %>%
        filter(Significance != "Not Significant") %>%
        top_n(10, -!!sym(pval_col)),  # p_val_adj top10
      aes(label = gene),
      size = 3,
      box.padding = 0.35,
      point.padding = 0.5,
      color = "black",
      segment.size = 0.2,
      segment.alpha = 0.4
    )
  
  ggsave(filename, width = 8, height = 6, dpi = 1200)

  print("zero pval genes:")
  print(zero_pval_genes)
}

for(i in unique(markers$cluster)){
    marker_subset <- filter(markers, cluster == i)
    plot_volcano(data = marker_subset, coef_col = coef_col, pval_col = pval_col, coef_threshold = coef_threshold, pval_threshold = pval_threshold, title = paste0("Volcano Plot in ", i), filename = paste0("volcano_plot_", i, ".pdf"))
}