# https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/computing_moments.ipynb
# This is only for development purposes

import sys
sys.path.append('/home/ssm-user/Github/scrna-parameter-estimation/dist/memento-0.0.5-py3.8.egg')
import memento
import matplotlib.pyplot as plt
import scanpy as sc
import memento
import pandas as pd
import seaborn as sns
import pickle as pkl

adata = sc.read(input_h5ad) # adata.X == raw counts
adata = adata[adata.obs.cell == 'CD14+ Monocytes'].copy()
print(adata)

adata.obs['stim'] = adata.obs['stim'].apply(lambda x: 0 if x == 'ctrl' else 1)
adata.obs[['ind', 'stim', 'cell']].sample(5)

# Create groups for hypothesis testing and compute 1D parameters
from scipy.sparse.csr import csr_matrix
type(adata.X) == csr_matrix

adata.obs['capture_rate'] = 0.07
memento.setup_memento(adata, q_column='capture_rate')

adata.obs.columns

memento.create_groups(adata, label_columns=['stim', 'ind'])
memento.compute_1d_moments(adata,
    min_perc_group=.9) # percentage of groups that satisfy the condition for a gene to be considered. 

# For a given gene identified, extract mean and residual variance estimates
mean, var, counts = memento.get_1d_moments(adata)
mean.query('gene == "IFI6"')
var.query('gene == "IFI6"')

# Create a figure with estimated moments
ctrl_cols = [c for c in mean.columns if '^0^' in c]
stim_cols = [c for c in mean.columns if '^1^' in c]

mean_ctrl, mean_stim = mean.query('gene == "IFI6"')[ctrl_cols].values.reshape(-1), mean.query('gene == "IFI6"')[stim_cols].values.reshape(-1)
var_ctrl, var_stim = var.query('gene == "IFI6"')[ctrl_cols].values.reshape(-1), var.query('gene == "IFI6"')[stim_cols].values.reshape(-1)

#make longform DataFrame for plotting. There are many ways to do this, this just lets you use seaborn in a straightforward way
df1 = pd.DataFrame()
df1['mean'] = mean_ctrl
df1['var'] = var_ctrl
df1['condition'] = 'ctrl'
df2 = pd.DataFrame()
df2['mean'] = mean_stim
df2['var'] = var_stim
df2['condition'] = 'stim'
df = pd.concat([df1, df2])

plt.figure(figsize=(5,3))
plt.subplots_adjust(wspace=0.5)
plt.subplot(1, 2, 1)
sns.boxplot(x='condition', y='mean',data=df)
sns.stripplot(x='condition', y='mean',data=df, edgecolor='gray', s=6, linewidth=1.5)
plt.subplot(1, 2, 2)
sns.boxplot(x='condition', y='var',data=df)
sns.stripplot(x='condition', y='var',data=df, edgecolor='gray', s=6, linewidth=1.5)
