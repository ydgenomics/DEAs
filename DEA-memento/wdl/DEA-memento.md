# Using [memento](https://github.com/yelabucsf/scrna-parameter-estimation) do DEA
- **Brief:** 一种考虑到分组间表达方差/多样性的差异分析软件，使用超几何检验对数据进行了一定修复，使其更加符合真实的数据分布
- **Fature** 优化可视化，提高结果的可解释性
- **Log:** 
  - v1.0.1 250813 035项目第一次优化

---
# Input
- **Variable**
  - 

---
# Output
- **Interpretation**
  Differential Variability Genes (DV) 的含义：差异变异性基因（Differential Variability Genes, DV）是指在不同条件下，基因表达的变异性存在显著差异的基因。具体来说，这些基因在某一条件下的表达变异性可能显著高于或低于其他条件下的变异性。这可能意味着这些基因在不同条件下的表达调控机制存在差异，或者受到不同因素的影响，导致其表达的稳定性不同。例如，在某种疾病状态下，某些基因的表达变异性可能增加，这可能与疾病的发生发展或异质性有关。

蓝色是de, 红色是dv
---
# Detail
- **Overview**

- **Software**
  - `memento` is a Python package for estimating the mean, variability, and gene correlation from scRNA-seq data as well as contructing a framework for hypothesis testing of differences in these parameters between groups of cells. Method-of-moments estimators are used for parameter estimation, and efficient resampling is used to construct confidence intervals and establish statistical significance.

- **Image**

---
# Reference & Citation
> [Cell || Resource || 2024 || Memento: 用于单细胞 RNA 测序数据差异表达分析的矩量框架方法](https://mp.weixin.qq.com/s/WGl51Y8DiIQHybq4cM_bRA)
> [差异基因找的不好？Cell刚发的这个单细胞差异统计的方法，可以用到咱们自己的数据上](https://mp.weixin.qq.com/s/5kwfiapfBKVSQmo5Zuh6EA)
> [【综述】Nature Methods | 干货！一文读懂单细胞转录组分析的现状和问题！](https://mp.weixin.qq.com/s/pF95r2p0KK1LSW-l5YFrsw)
> [单细胞转录组差异分析的8大痛点](https://mp.weixin.qq.com/s/4VThQEdclByz6jaJ3EQOdQ)
> [单细胞下游分析 | 基础知识 | ③差异分析](https://mp.weixin.qq.com/s/2sQ2cbEQv2VWmKpXCTmTrQ)
> [特异性基因筛选方法比较](https://mp.weixin.qq.com/s/8g4lSm8az7dyo__RZALL1w)


---
# Coder
- **Editor:** yangdong (yangdong@genomics.cn)
- **GitHub:** [ydgenomics](https://github.com/ydgenomics)
- **Prospect:** Focused on innovative, competitive, open-source projects and collaboration
- **Repository:** [DEAs/DEA-memento]()
---
