# DEA-Seurat
使用Seurat内置的三个函数做基因筛选/差异表达分析
输出的csv统一都包含这几列：`gene`,`cluster`,`p_val_adj`,`avg_log2FC`；便于后续的分析
[方法比较和结果解读：细胞类群marker基因识别及可视化](https://mp.weixin.qq.com/s/XA0gP-uYJmgcSQ1VAAYxYA)
[五种方式可视化Marker基因](https://mp.weixin.qq.com/s/iC8oB1LJD3y6y6sI-jdeQg)

[FindAllMarkers()](https://satijalab.org/seurat/reference/findallmarkers)
| p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | cluster | gene |
|-------|------------|-------|-------|-----------|---------|------|
[FindConservedMarkers()](https://satijalab.org/seurat/reference/findconservedmarkers)
| WT_p_val | WT_avg_log2FC | WT_pct.1 | WT_pct.2 | WT_p_val_adj | Mut_p_val | Mut_avg_log2FC | Mut_pct.1 | Mut_pct.2 | Mut_p_val_adj | max_pval | minimump_p_val | avg_log2FC | p_val_adj | cluster | gene |
|----------|---------------|----------|----------|--------------|-----------|----------------|-----------|-----------|---------------|----------|----------------|------------|-----------|---------|------|
[FindMarkers()](https://satijalab.org/seurat/reference/findmarkers)
| p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | cluster | gene_id | gene |
|-------|------------|-------|-------|-----------|---------|---------|------|

对于一些数据不仅包含两组对照数据，比如时序数据，可能会存在多个处理例如(T1,T2,T3,T4,T5,T6)，对于这种数据，我们要是关注某一个cluster中不同处理的差别话，可以对该cluster取subset，然后使用FindAllMarkers()，这个时候`group.by`就是不同时序，而不是cluster。