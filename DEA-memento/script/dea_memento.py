### Date; 250816 dea_memento.py
### Ref: https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/binary_testing_replicates.ipynb
### adata: stim(control) ind cell
### stim=biosample: s1, s2; ind=sample; cell=leiden_res_0.50: 0
### adata.X 一定要是原始数据且为CSR matrix

import memento
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import click
from scipy.sparse.csr import csr_matrix
import numpy as np
import os
import seaborn as sns
import itertools
import argparse
from adjustText import adjust_text

parser = argparse.ArgumentParser(description="DEA memento script parameters.")
parser.add_argument('--input_h5ad', type=str, default='/data/work/Single-Cell-Pipeline/DEA-memento/peanut.h5ad', help='Path to the input h5ad file.')
parser.add_argument('--group_key', type=str, default='biosample', help="Key for grouping samples.")
parser.add_argument('--control_value', type=str, default='H1314', help='Value of control group.')
parser.add_argument('--sample_key', type=str, default='sample', help='Key for sample identification.')
parser.add_argument('--capture_rate', type=float, default=0.07, help='Capture rate.')
parser.add_argument('--min_perc', type=float, default=0.7, help='Minimum percentage of groups that satisfy the condition for a gene to be considered.')
parser.add_argument('--pval_threshold', type=float, default=0.01, help='P-value threshold.')
parser.add_argument('--n_cpu', type=int, default=4, help='Number of CPUs to use for parallel processing.')
parser.add_argument('--top_number', type=int, default=10, help='Number of top genes to display.')
parser.add_argument('--perform_2d_test', type=str, default='no', help='Perform 2D hypothesis testing?')

args = parser.parse_args()

input_h5ad=args.input_h5ad
group_key=args.group_key
control_value=args.control_value
sample_key=args.sample_key
capture_rate=args.capture_rate
min_perc=args.min_perc
pval_threshold=args.pval_threshold
n_cpu=args.n_cpu
top_number=args.top_number
perform_2d_test=args.perform_2d_test

adata = sc.read(input_h5ad)
print(adata)
prefix = os.path.splitext(os.path.basename(input_h5ad))[0]
print(f"Prefix: {prefix}")

adata.obs[group_key] = adata.obs[group_key].apply(lambda x: 0 if x == control_value else 1)
adata.obs[[group_key, sample_key]].sample(5)
type(adata.X) == csr_matrix
if not isinstance(adata.X, csr_matrix):
    adata.X = csr_matrix(adata.X)
    

adata.obs['capture_rate'] = capture_rate
memento.setup_memento(adata, q_column='capture_rate')
memento.create_groups(adata, label_columns=[group_key, sample_key])
memento.compute_1d_moments(adata,
    min_perc_group=min_perc) # 0.7(default) percentage of groups that satisfy the condition for a gene to be considered. 

# ------------------------- Perform 1D hypothesis testing
sample_meta = memento.get_groups(adata)
sample_meta[sample_key] = sample_meta[sample_key].astype('category') # make sure to not confuse ourselves in case replicate labels are numbers
sample_meta.head(3)
treatment_df = sample_meta[[group_key]]
treatment_df.head(5)
cov_df = pd.get_dummies(sample_meta[sample_key].astype('category'))
cov_df.head(3)
memento.ht_1d_moments(
    adata, 
    treatment=treatment_df,
    covariate=cov_df,
    num_boot=5000, 
    verbose=1,
    num_cpus=n_cpu)
result_1d = memento.get_1d_ht_result(adata)
result_1d['type'] = "N"

# -------------------- plot -----------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text   # 新增

# -------------------------------------------------------------

# 1) 标记 Signature
sig_mask = (result_1d['dv_pval'] < pval_threshold) | (result_1d['de_pval'] < pval_threshold)
result_1d.loc[sig_mask, 'type'] = 'Signature'
colors = np.where(result_1d['type'] == 'Signature', 'black', 'gray')

# 2) 要标注的基因
top_dv = result_1d.nsmallest(top_number, 'dv_pval')
top_de = result_1d.nsmallest(top_number, 'de_pval')
label_df = pd.concat([top_dv, top_de]).drop_duplicates(subset='gene')

# -------------------------------------------------------------
# 3) 画图
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(result_1d['de_coef'], result_1d['dv_coef'],
           s=5, c=colors, alpha=0.8, zorder=2)

texts = []
# dv_pval Top5：红色
for _, row in top_dv.iterrows():
    texts.append(
        ax.text(row['de_coef'], row['dv_coef'], row['gene'],
                fontsize=6, color='red')
    )

# de_pval Top5：蓝色
for _, row in top_de.iterrows():
    texts.append(
        ax.text(row['de_coef'], row['dv_coef'], row['gene'],
                fontsize=6, color='blue')
    )

# 自动避免重叠 + 画指引虚线
adjust_text(texts,
            arrowprops=dict(arrowstyle='-', color='black', lw=0.6, ls='--'),
            ax=ax)
ax.axhline(0, color='gray', lw=0.8, zorder=1)   # 水平线 y=0
ax.axvline(0, color='gray', lw=0.8, zorder=1)   # 垂直线 x=0
ax.set_title(f'Differential analysis for {prefix}')
ax.set_xlabel('Differential Expression Coefficient')
ax.set_ylabel('Differential Variability Coefficient')
plt.tight_layout()
plt.savefig(f'differential_expression_replicate_{prefix}.pdf')
plt.close()
result_1d.to_csv(f'result_1d_replicate_{prefix}.txt', sep='\t', index=False)

def plot_top_genes(cell_adata, top_genes, cell_type_dir, perform_2d_test, cell_type, group_key, sample_key):
    os.makedirs('topgenes_boxplot', exist_ok=True)
    mean, var, counts = memento.get_1d_moments(cell_adata)
    for gene in top_genes:
        plot_gene_moments(mean, var, gene, cell_type_dir)
        os.chdir('..')
    if perform_2d_test == 'yes':
        os.makedirs('topgenes_2d_test', exist_ok=True)
        for gene in top_genes:
            perform_2d_hypothesis_testing(cell_adata, gene, cell_type_dir, cell_type, group_key, sample_key)
            os.chdir('..')

def plot_gene_moments(mean, var, gene, cell_type_dir):
    mean.query(f'gene == "{gene}"')
    var.query(f'gene == "{gene}"')
    ctrl_cols = [c for c in mean.columns if '^0^' in c]
    stim_cols = [c for c in mean.columns if '^1^' in c]
    mean_ctrl, mean_stim = mean.query(f"gene == '{gene}'")[ctrl_cols].values.reshape(-1), mean.query(f"gene == '{gene}'")[stim_cols].values.reshape(-1)
    var_ctrl, var_stim = var.query(f"gene == '{gene}'")[ctrl_cols].values.reshape(-1), var.query(f"gene == '{gene}'")[stim_cols].values.reshape(-1)
    os.chdir(os.path.join(cell_type_dir, 'topgenes_boxplot'))
    #make longform DataFrame for plotting. There are many ways to do this, this just lets you use seaborn in a straightforward way
    df1 = pd.DataFrame({'mean': mean_ctrl, 'var': var_ctrl, 'condition': 'ctrl'})
    df2 = pd.DataFrame({'mean': mean_stim, 'var': var_stim, 'condition': 'stim'})
    df = pd.concat([df1, df2])
    plt.figure(figsize=(6,4))
    plt.subplots_adjust(wspace=0.5)
    plt.subplot(1, 2, 1)
    sns.boxplot(x='condition', y='mean', data=df)
    sns.stripplot(x='condition', y='mean', data=df, edgecolor='gray', s=6, linewidth=1.5)
    plt.subplot(1, 2, 2)
    sns.boxplot(x='condition', y='var', data=df)
    sns.stripplot(x='condition', y='var', data=df, edgecolor='gray', s=6, linewidth=1.5)
    plt.savefig(f'output_{gene}_boxplot.pdf')
    plt.close()

def perform_2d_hypothesis_testing(cell_adata, gene, cell_type_dir, cell_type, group_key, sample_key):
    os.chdir(os.path.join(cell_type_dir, 'topgenes_2d_test'))
    sample_meta = memento.get_groups(cell_adata)
    sample_meta[sample_key] = sample_meta[sample_key].astype('category')
    treatment_df = sample_meta[[group_key]]
    cov_df = pd.get_dummies(sample_meta[sample_key].astype('category'))
    #result_2d.to_csv(f'result_2d_replicate_{cell_adata.obs["cell"].unique()[0]}_{gene}.txt', sep='\t', index=False)
    gene_pairs = list(itertools.product([gene], cell_adata.var.index.tolist()))
    memento.compute_2d_moments(cell_adata, gene_pairs)
    memento.ht_2d_moments(cell_adata, treatment=treatment_df, covariate=cov_df, num_boot=5000, verbose=1, num_cpus=n_cpu)
    result_2d = memento.get_2d_ht_result(cell_adata)
    result_2d.to_csv(f'result_2d_replicate_{cell_type}_{gene}.txt', sep='\t', index=False)

gene_de = result_1d.sort_values(by='de_pval').head(top_number)['gene']
gene_dv = result_1d.sort_values(by='dv_pval').head(top_number)['gene']
plot_top_genes(adata, gene_de, cell_type_dir= "./", perform_2d_test=perform_2d_test, cell_type=prefix, group_key=group_key, sample_key=sample_key)
plot_top_genes(adata, gene_dv, cell_type_dir= "./", perform_2d_test=perform_2d_test, cell_type=prefix, group_key=group_key, sample_key=sample_key)
