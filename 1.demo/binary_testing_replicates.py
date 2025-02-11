# https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/binary_testing_replicates.ipynb
# This is only for development purposes

import sys
import memento
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

adata = sc.read(input_h5ad) # adata.X == raw counts
adata = adata[adata.obs.cell == 'CD14+ Monocytes'].copy()
print(adata)

adata.obs['stim'] = adata.obs['stim'].apply(lambda x: 0 if x == 'ctrl' else 1)
adata.obs[['ind', 'stim', 'cell']].sample(5)

# Details on designating replicates
from scipy.sparse.csr import csr_matrix

type(adata.X) == csr_matrix

adata.obs['capture_rate'] = 0.07
memento.setup_memento(adata, q_column='capture_rate')

memento.create_groups(adata, label_columns=['stim', 'ind'])

memento.compute_1d_moments(adata,
    min_perc_group=.9) # percentage of groups that satisfy the condition for a gene to be considered. 

# Perform 1D hypothesis testing
sample_meta = memento.get_groups(adata)
sample_meta['ind'] = sample_meta['ind'].astype('category') # make sure to not confuse ourselves in case replicate labels are numbers.
sample_meta.head(3)

treatment_df = sample_meta[['stim']]
treatment_df.head(5)

cov_df = pd.get_dummies(sample_meta['ind'].astype('category'))
cov_df.head(3)

memento.ht_1d_moments(
    adata, 
    treatment=treatment_df,
    covariate=cov_df,
    num_boot=5000, 
    verbose=1,
    num_cpus=40)

result_1d = memento.get_1d_ht_result(adata)

import numpy as np

plt.scatter(result_1d.de_coef, result_1d.dv_coef, s=1)

import scipy.stats as stats

result_1d.query('de_coef > 0').sort_values('de_pval')
result_1d.query('dv_coef > 0 & de_coef > 0').sort_values('dv_pval').head(10)

# Perform 2D hypothesis testing
import itertools

gene_pairs = list(itertools.product(['IRF7'], adata.var.index.tolist()))

memento.compute_2d_moments(adata, gene_pairs)

memento.ht_2d_moments(
    adata, 
    treatment=treatment_df,
    covariate=cov_df,
    num_boot=5000, 
    verbose=1,
    num_cpus=40)

result_2d = memento.get_2d_ht_result(adata)
result_2d
result_2d.sort_values('corr_pval').head(10)