# https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/binary_testing.ipynb
# input_key: input_h5ad, input_key,
import memento
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

adata = sc.read(input_h5ad) # adata.X == raw counts
adata = adata[adata.obs['cell'] == 'CD14+ Monocytes'].copy()
print(adata)

adata.obs['stim'] = adata.obs['stim'].apply(lambda x: 0 if x == 'ctrl' else 1)
adata.obs[['ind', 'stim', 'cell']].sample(5)

# Perform the binary comparison between two groups
result_1d = memento.binary_test_1d(
    adata=adata, 
    capture_rate=0.07, 
    treatment_col='stim', 
    num_cpus=12,
    num_boot=5000)

plt.scatter(result_1d.de_coef, result_1d.dv_coef, s=1)

# Print the top differential mean genes (DM) and top differential variability genes (DV)
result_1d.query('de_coef > 0').sort_values('de_pval').head(10)
result_1d.query('dv_coef > 0 & de_coef > 0').sort_values('dv_pval').head(10)

# Perform 2D hypothesis testing
import itertools

gene_pairs = list(itertools.product(['IRF7'], adata.var.index.tolist()))

result_2d = memento.binary_test_2d(
    adata=adata, 
    gene_pairs=gene_pairs, 
    capture_rate=0.07, 
    treatment_col='stim', 
    num_cpus=12, 
    num_boot=5000)

result_2d.sort_values('corr_pval').head(5)