# DEA-Seurat
使用Seurat内置的三个函数做基因筛选/差异表达分析
  - FindAllMarkers适用于分群后找各个群的marker基因，该群区别于其它群特异的基因(pos/neg) `assay` `group.by` `only.pos = TRUE`
  - FindConservedMarkers适用于整合后分群找各个群的marker基因，该群应该由多个批次数据组成故要找到能够代表多个批次的该群marker基因 `assay` `ident.1` `grouping.var` `only.pos = TRUE`
  - Findmarkers使用于两分组找差异基因 `assay` `ident.1` `ident.2` `group.by`

# Q&A
1. 找差异基因应该包括差异高表达基因和差异低表达基因，而差异低表达基因很可能受限于测序精度，真实的生物学上可能是有表达的，所以我们更加关注于差异高表达基因，更加的可靠，应设置`only.pos = TRUE`，默认参数是`FALSE`。其中一个问题就是在绘制火山图时，选了only.pos那么自然低表达基因会少很多。
2. 面向的数据是标准后还是原始数据呢？memento是以counts原始数据作为输入，而Seurat的这几个参数不太清楚。个人认为如果单一实验内cluster相比，用counts和data的效果应该都是一样的；如果是两个实验集，如果捕获差距不大，也可以认为是同一实验；但是如果捕获差距大，尤其是来自不同技术的数据，可能用data要好一点？

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

对于一些数据不仅包含两组对照数据，比如时序数据，可能会存在多个处理例如(T1,T2,T3,T4,T5,T6)，对于这种数据，我们要是关注某一个cluster中不同处理的差别话，可以对该cluster取subset，然后使用FindAllMarkers()，这个时候`group.by`就是不同时序，而不是cluster，找到的就应该是该细胞类型的时期特异基因。