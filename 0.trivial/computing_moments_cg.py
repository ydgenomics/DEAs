# https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/computing_moments.ipynb
# This is only for development purposes

import sys
import memento
import matplotlib.pyplot as plt
import scanpy as sc
import memento
import pandas as pd
import seaborn as sns
import pickle as pkl
from scipy.sparse.csr import csr_matrix
import click

@click.command()
@click.option('--input_h5ad', default='/data/work/input/memento/Pog_unsoupx_memento.h5ad', help='Path to the input h5ad file.')
@click.option('--ctrl_name', default='V2.5R2404290045', help='Control sample name.')
@click.option('--cell_type_list', default='0', help='Comma-separated list of cell types.')
@click.option('--research_genes', default='P.cirratum_15658.t1', help='Comma-separated list of research genes.')
@click.option('--treatment_col', default='stim', help='Column name for treatment.')


def main(input_h5ad, ctrl_name, cell_type_list, research_genes, treatment_col):
    cell_type_list = cell_type_list.split(',')
    research_genes = research_genes.split(',')

    adata = sc.read(input_h5ad) # adata.X == raw counts
    #adata.X = adata.layers["counts"].copy() # adata.X == raw counts
    print(adata)
    # Loop through unique cell types in adata.obs['cell']
    for cell_type in cell_type_list:
        print(f"Processing cell type: {cell_type}")
        
        # Filter the data for the current cell type
        cell_adata = adata[adata.obs['cell'] == cell_type].copy()
        cell_adata.obs['stim'] = cell_adata.obs['stim'].apply(lambda x: 0 if x == ctrl_name else 1)
        cell_adata.obs[['ind', 'stim', 'cell']].sample(5)
        print(cell_adata.obs['cell'].unique())

        # Create groups for hypothesis testing and compute 1D parameters
        type(cell_adata.X) == csr_matrix
        cell_adata.obs['capture_rate'] = 0.07
        memento.setup_memento(cell_adata, q_column='capture_rate')
        cell_adata.obs.columns

        memento.create_groups(cell_adata, label_columns=['stim', 'ind'])
        memento.compute_1d_moments(cell_adata,
                                   min_perc_group=.9) # percentage of groups that satisfy the condition for a gene to be considered. 

        # For a given gene identified, extract mean and residual variance estimates
        mean, var, counts = memento.get_1d_moments(cell_adata)
        for gene in research_genes:
            mean.query(f"gene == '{gene}'")
            var.query(f"gene == '{gene}'")
            # Create a figure with estimated moments
            ctrl_cols = [c for c in mean.columns if '^0^' in c]
            stim_cols = [c for c in mean.columns if '^1^' in c]
            
            mean_ctrl, mean_stim = mean.query(f"gene == '{gene}'")[ctrl_cols].values.reshape(-1), mean.query(f"gene == '{gene}'")[stim_cols].values.reshape(-1)
            var_ctrl, var_stim = var.query(f"gene == '{gene}'")[ctrl_cols].values.reshape(-1), var.query(f"gene == '{gene}'")[stim_cols].values.reshape(-1)

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
            plt.savefig(f'output_{gene}_boxplot.pdf')

if __name__ == '__main__':
    main()
