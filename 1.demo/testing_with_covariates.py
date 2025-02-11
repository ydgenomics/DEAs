# https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/testing_with_covariates.ipynb
# This is only for development purposes

import sys
# sys.path.append('/home/ssm-user/Github/scrna-parameter-estimation/dist/memento-0.0.8-py3.8.egg')
import memento
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import pickle as pkl

adata = sc.read(input_h5ad) # adata.X == raw counts
adata = adata[adata.obs.cell == 'CD14+ Monocytes'].copy()
print(adata)

adata.obs['stim'] = adata.obs['stim'].apply(lambda x: 0 if x == 'ctrl' else 1)
adata.obs[['ind', 'stim', 'cell']].sample(5)

# Engineer the covariate to be used in memento.
# These are not actually chromosome 1 genes
# TODO: Remake this tutorial with actual chr1 labels and maybe with the aneuploidy dataset so that it makes sense
chr1_genes = list(np.random.choice(adata.var.index, 4000))

adata_chrom = adata.copy().copy()
adata_chrom.obs['chr_expr'] = adata_chrom[:, chr1_genes].X.sum(axis=1).astype(int)
adata_chrom.obs['chr_expr_bin'] = pd.qcut(adata_chrom.obs['chr_expr'], 10)
adata_chrom.obs = adata_chrom.obs.join(adata_chrom.obs.groupby('chr_expr_bin')['chr_expr'].median(), on='chr_expr_bin', rsuffix='_avg')

adata_chrom.obs.head(2)

# Setup memento with the treatment and covariates
adata_chrom.obs['capture_rate'] = 0.07
memento.setup_memento(adata_chrom, q_column='capture_rate')
memento.create_groups(adata_chrom, label_columns=['stim', 'chr_expr_avg'])
memento.compute_1d_moments(adata_chrom, min_perc_group=.7, gene_list=chr1_genes)

adata_chrom.shape

# Perform 1D hypothesis testing
sample_meta = memento.get_groups(adata_chrom)

# The covariate DataFrame - pick the covariate columns
cov_df = sample_meta[['chr_expr_avg']]

# The treatment DataFrame - pick the treatment column
treat_df = sample_meta[['stim']]

memento.ht_1d_moments(
    adata_chrom, 
    treatment=treat_df,
    covariate=cov_df,
    resampling='bootstrap',
    num_boot=5000, 
    verbose=1,
    num_cpus=14)

result_1d = memento.get_1d_ht_result(adata_chrom)
result_1d.query('de_coef > 0').sort_values('de_pval').head(10)


