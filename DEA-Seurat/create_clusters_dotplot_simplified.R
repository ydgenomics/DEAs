# Title: create_clusters_dotplot_simplified.R
# Author: husasa, ydgenomics
# Date: 2025-05-14
# Csv need these columns: gene, cluster, p_val_adj, avg_log2FC

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(optparse)

option_list <- list(
  make_option(c("-g", "--gene_csv"), type = "character", default = "/data/work/test/result/markers_peanut.csv",
              help = "Input file path, default is '/data/work/test/result/markers_peanut.csv'")
)
opt <- parse_args(OptionParser(option_list = option_list))
markers <- read.csv(opt$gene_csv, row.names = NULL)
head(markers)
markers$cell_types <- markers$cluster
cell_types <- unique(markers$cell_types)
# 添加显著性标签
markers$label <- ifelse(markers$p_val_adj < 0.01, "adjust P-val<0.01", "adjust P-val>=0.01")
head(markers)

# 添加调控方向
markers$significant <- "no"
markers$significant[markers$p_val_adj < 0.01 & markers$avg_log2FC > 0.5] <- "sigUp"
markers$significant[markers$p_val_adj < 0.01 & markers$avg_log2FC < -0.5] <- "sigDown"

# 合并数据
alldiff_clusters <- markers

# 计算并打印显著差异基因数量
sig_up <- sum(markers$significant == "sigUp")
sig_down <- sum(markers$significant == "sigDown")
cat(paste0("显著上调基因数量: ", sig_up, "\n"))
cat(paste0("显著下调基因数量: ", sig_down, "\n"))



# 将细胞类型设为因子变量，保持指定顺序
alldiff_clusters$cell_type <- factor(alldiff_clusters$cell_type, levels = cell_types)

# 提取显著差异基因
sig_genes <- alldiff_clusters %>% 
  filter(significant != "no")

# 获取每个细胞类型中avg_log2FC的最大值和最小值，用来画背景
# 获取最大值
celltype_avg_log2FC_max <- c()
for (i in 1:length(cell_types)) {
  x <- cell_types[i]
  current_data <- filter(alldiff_clusters, cell_type == x, significant != "no")
  if (nrow(current_data) > 0) {
    current_max <- max(current_data$avg_log2FC) 
  } else {
    current_max <- 0
  }
  celltype_avg_log2FC_max <- append(celltype_avg_log2FC_max, current_max)
}
celltype_avg_log2FC_max <- as.numeric(celltype_avg_log2FC_max)

# 获取最小值
celltype_avg_log2FC_min <- c()
for (i in 1:length(cell_types)) {
  x <- cell_types[i]
  current_data <- filter(alldiff_clusters, cell_type == x, significant != "no")
  if (nrow(current_data) > 0) {
    current_min <- min(current_data$avg_log2FC) 
  } else {
    current_min <- 0
  }
  celltype_avg_log2FC_min <- append(celltype_avg_log2FC_min, current_min)
}
celltype_avg_log2FC_min <- as.numeric(celltype_avg_log2FC_min)

# 将最大值和最小值构建数据框
dfbar <- data.frame(
  x = cell_types,
  y = celltype_avg_log2FC_max
)
dfbar$x <- factor(dfbar$x, levels = cell_types)

dfbar1 <- data.frame(
  x = cell_types,
  y = celltype_avg_log2FC_min
)
dfbar1$x <- factor(dfbar1$x, levels = cell_types)

# 创建显示上调和下调基因数量的统计数据
sig_counts <- alldiff_clusters %>%
  group_by(cell_type, significant) %>%
  filter(significant != "no") %>%
  summarize(count = n(), .groups = 'drop')

# 计算每个细胞类型对应的标签和基因数量
dfcol <- data.frame(
  x = cell_types,
  y = 0,
  label = cell_types,  # 使用细胞类型名称作为标签
  count_up = sapply(cell_types, function(ct) {
    sum(sig_counts$count[sig_counts$cell_type == ct & sig_counts$significant == "sigUp"])
  }),
  count_down = sapply(cell_types, function(ct) {
    sum(sig_counts$count[sig_counts$cell_type == ct & sig_counts$significant == "sigDown"])
  })
)
dfcol$x <- factor(dfcol$x, levels = cell_types)

# 为每个细胞类型选择不同的颜色
n_clusters <- length(cell_types)
cell_type_colors <- scales::hue_pal()(n_clusters)

# 使用用户指定的颜色
up_color <- "#F8766D"    # 上调基因颜色
down_color <- "#00BFC4"  # 下调基因颜色

# 获取每个分组中表达差异最显著的Top2上调基因和Top2下调基因
top2_up_list <- vector("list", length(cell_types))
top2_down_list <- vector("list", length(cell_types))

for (i in 1:length(cell_types)) {
  x <- cell_types[i]
  # 获取Top2上调基因
  current_top_up <- filter(alldiff_clusters, cell_type == x, significant == "sigUp") %>%
    distinct(gene, .keep_all = TRUE) %>% 
    top_n(2, abs(avg_log2FC))  # 选择Top2上调基因
  top2_up_list[[i]] <- current_top_up
  
  # 获取Top2下调基因
  current_top_down <- filter(alldiff_clusters, cell_type == x, significant == "sigDown") %>% 
    distinct(gene, .keep_all = TRUE) %>% 
    top_n(2, abs(avg_log2FC))  # 选择Top2下调基因
  top2_down_list[[i]] <- current_top_down
}

top2_up_all <- do.call(rbind, top2_up_list)
top2_down_all <- do.call(rbind, top2_down_list)

# 合并上调和下调的Top2基因
top2_all <- rbind(top2_up_all, top2_down_all)

# 将细胞类型设为因子变量，保持指定顺序
if (nrow(top2_all) > 0) {
  top2_all$cell_type <- factor(top2_all$cell_type, levels = cell_types)
}

# 绘制图形
p <- ggplot() +
  # 灰色背景柱 - 上调
  geom_col(data = dfbar,
           mapping = aes(x = x, y = y),
           fill = "#f0f0f0", alpha = 0.6) +
  # 灰色背景柱 - 下调
  geom_col(data = dfbar1,
           mapping = aes(x = x, y = y),
           fill = "#f0f0f0", alpha = 0.6) +
  # 在底部添加彩色标记，高度为1.3
  geom_tile(data = dfcol,
            aes(x = x, y = y, fill = x),
            height = 1.3,  # 高度设为1.3
            show.legend = FALSE) +
  scale_fill_manual(values = cell_type_colors) +
  # 显著差异基因散点 - 上调基因（向上移动点位置）
  geom_jitter(data = sig_genes %>% filter(significant == "sigUp"),
              aes(x = cell_type, y = avg_log2FC, color = significant),
              size = 1.5,
              position = position_jitter(width = 0.4, height = 0.01, seed = 123)) +
  # 显著差异基因散点 - 下调基因（向下移动点位置）
  geom_jitter(data = sig_genes %>% filter(significant == "sigDown"),
              aes(x = cell_type, y = avg_log2FC, color = significant),
              size = 1.5,
              position = position_jitter(width = 0.4, height = 0.01, seed = 456)) +
  # Top2显著差异基因 - 上调（点更大）
  geom_jitter(data = top2_up_all,
              aes(x = cell_type, y = avg_log2FC, color = significant),
              size = 4,
              position = position_jitter(width = 0.4, height = 0.01, seed = 789)) +
  # Top2显著差异基因 - 下调（点更大）
  geom_jitter(data = top2_down_all,
              aes(x = cell_type, y = avg_log2FC, color = significant),
              size = 4,
              position = position_jitter(width = 0.4, height = 0.01, seed = 987)) +              
  # 使用指定的自定义颜色
  scale_color_manual(values = c("sigUp" = up_color, "sigDown" = down_color, "no" = "grey80"),
                     labels = c("sigUp" = "上调", "sigDown" = "下调", "no" = "无显著差异")) +
  # 给Top2基因添加文本标签
  geom_text_repel(
    data = top2_all,
    aes(x = cell_type, y = avg_log2FC, label = gene),
    force = 1.5,  # 增加力度，让标签更分散
    box.padding = 0.7,  # 增加空间
    point.padding = 0.5,
    max.overlaps = 20,
    size = 3.5,
    arrow = arrow(length = unit(0.01, "inches"), type = "open", ends = "last"),
    family = "sans"
  ) +    
  # 在每个细胞类型上方显示上调基因数量
  geom_text(
    data = dfcol,
    aes(x = x, y = max(celltype_avg_log2FC_max) + 0.5, label = count_up),
    color = up_color,
    size = 3.5,
    family = "sans"
  ) +
  # 在每个细胞类型下方显示下调基因数量
  geom_text(
    data = dfcol,
    aes(x = x, y = min(celltype_avg_log2FC_min) - 0.5, label = count_down),
    color = down_color,
    size = 3.5,
    family = "sans"
  ) +
  # 设置图表标题和轴标签
  labs(
    title = "Clusters",
    x = "",
    y = "Average log2 fold change"
  ) +
  # 美化主题
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "sans"),
    axis.title.y = element_text(size = 14, family = "sans"),
    axis.text.y = element_text(size = 12, family = "sans"),
    axis.text.x = element_blank(),
    axis.line.y = element_line(color = "black", size = 1),
    axis.line.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, family = "sans")
  )

# 保存图形
ggsave("simplified_clusters_differential_genes.pdf", p, width = 12, height = 8)
ggsave("simplified_clusters_differential_genes.png", p, width = 12, height = 8, dpi = 300)

cat("图形已保存为 simplified_clusters_differential_genes.pdf 和 simplified_clusters_differential_genes.png\n")

# 创建汇总表格
cat("\n========== 创建差异基因汇总表格 ==========\n")
summary_table <- data.frame(
  cluster = cell_types,
  total_genes = sapply(cell_types, function(ct) sum(alldiff_clusters$cell_type == ct)),
  sig_genes = sapply(cell_types, function(ct) sum(alldiff_clusters$cell_type == ct & alldiff_clusters$significant != "no")),
  up_genes = sapply(cell_types, function(ct) sum(alldiff_clusters$cell_type == ct & alldiff_clusters$significant == "sigUp")),
  down_genes = sapply(cell_types, function(ct) sum(alldiff_clusters$cell_type == ct & alldiff_clusters$significant == "sigDown"))
)

# 添加上调下调基因差值列
summary_table$diff_up_down <- summary_table$up_genes - summary_table$down_genes

# 保存汇总表格
write.csv(summary_table, "simplified_differential_genes_summary.csv", row.names = FALSE)
cat("差异基因汇总表格已保存为 simplified_differential_genes_summary.csv\n")

# 保存显著差异基因数据
write.csv(sig_genes, "simplified_significant_differential_genes.csv", row.names = FALSE)
cat("显著差异基因数据已保存为 simplified_significant_differential_genes.csv\n") 