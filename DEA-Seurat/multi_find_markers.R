# Title: multi_find_markers.R
# Date: 2025-05-16
# Image: plantphone-R-02
# https://github.com/ydgenomics/git_test/blob/main/SCRIPTs/multi_find_markers.R
# This script is designed to analyze single-cell RNA sequencing data using the Seurat package.
# It provides functionality to identify marker genes across different groups or conditions.
# The script supports multiple tools such as FindAllMarkers, FindConservedMarkers, FindMarkers, and FindMultiMarkers.
# Users can specify various parameters like assay type, slot, grouping variables, and thresholds via command-line options.
# The results are saved as CSV files for further analysis.
# Same and Important columns of Different methods : gene, cluster, p_val_adj, avg_log2FC


library(Seurat)
library(dplyr)
library(optparse)

option_list <- list(
    make_option(c("-r", "--rds"), type = "character", default = "/data/work/test/peanut_dataget_Anno_concat.cg.rds", help = "Path to RDS file"),
    make_option(c("-t", "--tool"), type = "character", default = "FindAllMarkers,FindConservedMarkers,FindMarkers,FindMultiMarkers", help = "Tools to use"),
    make_option(c("-a", "--assay"), type = "character", default = "RNA", help = "Assay to use"),
    make_option(c("-s", "--slot"), type = "character", default = "data", help = "Slot to use"),
    make_option(c("-b", "--batch_var"), type = "character", default = "biosample", help = "Batch variable"),
    make_option(c("-g", "--group_var"), type = "character", default = "cell", help = "Group variable"),
    make_option(c("-p", "--min_pct"), type = "numeric", default = 0.01, help = "Minimum percentage"),
    make_option(c("-l", "--log_fc"), type = "numeric", default = 0.1, help = "Log fold change threshold"),
    make_option(c("-n", "--name"), type = "character", default = "peanut", help = "Output name"),
    make_option(c("-1", "--batch_1"), type = "character", default = "WT", help = "Batch 1"),
    make_option(c("-2", "--batch_2"), type = "character", default = "Mut", help = "Batch 2")
)
opt <- parse_args(OptionParser(option_list = option_list))
rds <- opt$rds
tool <- opt$tool
assay <- opt$assay
slot <- opt$slot
sample_var <- opt$batch_var
group_var <- opt$group_var
min_pct <- opt$min_pct
log_fc <- opt$log_fc
name <- opt$name
sample_1 <- opt$batch_1
sample_2 <- opt$batch_2

seu <- readRDS(rds)
seu
tool <- strsplit(tool, ",")[[1]]
DefaultAssay(seu) <- assay
celltypes <- unique(seu@meta.data[[group_var]])

if (assay == "RNA") {
    seu <- NormalizeData(seu)
}

# FindAllMarkers 
# "p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj" "cluster" "gene"
if ("FindAllMarkers" %in% tool) {
    Idents(seu) <- group_var
    allmarkers <- FindAllMarkers(seu, assay = assay, slot = slot, group.by = group_var, only.pos = FALSE, min.pct = min_pct, logfc.threshold = log_fc) # nolint
    print(dim(allmarkers))
    print(head(allmarkers))
    write.csv(allmarkers, paste0("allmarkers_", name, ".csv"), row.names = TRUE)
}

# FindConservedMarkers
# "WT_p_val" "WT_avg_log2FC" "WT_pct.1" "WT_pct.2" "WT_p_val_adj" "Mut_p_val" "Mut_avg_log2FC" "Mut_pct.1" "Mut_pct.2" "Mut_p_val_adj" "max_pval" "minimump_p_val" "cluster" "gene" "avg_log2FC" "p_val_adj"
if ("FindConservedMarkers" %in% tool) {
    conserved_markers <- list()
    for (cell_type in celltypes) {
        Idents(seu) <- group_var
        markers <- FindConservedMarkers(
            seu, 
            ident.1 = cell_type, 
            grouping.var = sample_var, 
            assay = assay, 
            slot = slot
        ) # nolint
        avg_log2FC_columns <- grep("_avg_log2FC$", names(markers), value = TRUE)
        if (length(avg_log2FC_columns) == 1) {
            markers$avg_log2FC <- markers[[avg_log2FC_columns]]
        } else {
            markers$avg_log2FC <- rowMeans(markers[, avg_log2FC_columns], na.rm = TRUE)
        }
        p_val_adj_columns <- grep("_p_val_adj$", names(markers), value = TRUE)
        if (length(avg_log2FC_columns) == 1) {
            markers$p_val_adj <- markers[[p_val_adj_columns]]
        } else {
            markers$p_val_adj <- rowMeans(markers[, p_val_adj_columns], na.rm = TRUE)
        }
        markers$cluster <- cell_type
        conserved_markers[[cell_type]] <- markers
    }
    conserved_markers <- bind_rows(conserved_markers)
    conserved_markers$gene <- rownames(conserved_markers)
    print(head(conserved_markers))
    write.csv(conserved_markers, paste0("conserved_markers_", name, ".csv"), row.names = TRUE) # nolint
}

new_col_name <- paste0(sample_var, "_", group_var)
new_col_values <- paste(seu@meta.data[[sample_var]], seu@meta.data[[group_var]], sep = "_")
seu <- AddMetaData(seu, metadata = new_col_values, col.name = new_col_name)
head(seu@meta.data)

# FindMarkers
# "p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj" "cluster" "gene"
if ("FindMarkers" %in% tool) {
  markers <- list()
  for (cell_type in celltypes) {
    seu_subset <- subset(seu, subset = !!sym(group_var) == cell_type)
    if (length(unique(seu_subset@meta.data[[new_col_name]])) == length(unique(seu@meta.data[[sample_var]]))) {
        Idents(seu_subset) <- seu_subset@meta.data[[sample_var]]
        markers_subset <- FindMarkers(
          seu_subset,
          ident.1 = sample_1,
          ident.2 = sample_2
        )
        markers_subset$cluster <- cell_type
        markers[[cell_type]] <- markers_subset
    } else {
        print(paste0("Attention: At least one sample subset has zero rows.", cell_type))
    }
  }
  markers <- bind_rows(markers)
  markers$gene <- rownames(markers)
  print(head(markers))
  write.csv(markers,paste0("markers_", name, ".csv"),row.names = FALSE) # nolint
}

# FindMultiMarkers
# "p_val" "avg_log2FC" "pct.1" "pct.2" "p_val_adj" "cluster" "gene" "compare_cluster" 
if ("FindMultiMarkers" %in% tool) {
  multi_markers <- list()
  for (cell_type in celltypes) {
    seu_subset <- subset(seu, subset = !!sym(group_var) == cell_type)
    if (length(unique(seu_subset@meta.data[[new_col_name]])) == length(unique(seu@meta.data[[sample_var]]))) {
        Idents(seu_subset) <- seu_subset@meta.data[[sample_var]]
        markers_subset <- FindAllMarkers(
          seu_subset,
          assay = assay,
          slot = slot,
          group.by = sample_var,
          only.pos = FALSE,
          min.pct = min_pct,
          logfc.threshold = log_fc
        ) # nolint
        markers_subset$compare_cluster <- markers_subset$cluster
        markers_subset$cluster <- cell_type
        multi_markers[[cell_type]] <- markers_subset
    } else {
        print(paste0("Attention: At least one sample subset has zero rows.", cell_type))
    }
  }
  multi_markers <- bind_rows(multi_markers)
  print(head(multi_markers))
  write.csv(multi_markers, paste0("multi_markers_", name, ".csv"), row.names = TRUE)
}