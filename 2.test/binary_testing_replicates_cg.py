# https://nbviewer.org/github/yelabucsf/scrna-parameter-estimation/blob/master/tutorials/binary_testing_replicates.ipynb
# This is only for development purposes

import sys
import memento
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import click
from scipy.sparse.csr import csr_matrix
import numpy as np
import scipy.stats as stats
import itertools

@click.command()
@click.option('--input_h5ad', default='/data/work/input/memento/Pog_unsoupx_memento.h5ad', help='Path to the input h5ad file.')
@click.option('--ctrl_name', default='V2.5R2404290045', help='Control sample name.')
@click.option('--cell_type_list', default='0', help='Comma-separated list of cell types.')
#@click.option('--research_genes', default='P.cirratum_15658.t1', help='Comma-separated list of research genes.')
@click.option('--treatment_col', default='stim', help='Column name for treatment.')

def main(input_h5ad, ctrl_name, cell_type_list, treatment_col):
    cell_type_list = cell_type_list.split(',')
    research_genes = research_genes.split(',')

    adata = sc.read(input_h5ad) # adata.X == raw counts!
    print(adata)

    # Loop through unique cell types in adata.obs['cell']
    for cell_type in cell_type_list:
        print(f"Processing cell type: {cell_type}")
        
        # Filter the data for the current cell type
        cell_adata = adata[adata.obs['cell'] == cell_type].copy()
        cell_adata.obs['stim'] = cell_adata.obs['stim'].apply(lambda x: 0 if x == ctrl_name else 1)
        cell_adata.obs[['ind', 'stim', 'cell']].sample(5)
        print(cell_adata.obs['cell'].unique())
        
        # Details on designating replicates
        type(cell_adata.X) == csr_matrix

        cell_adata.obs['capture_rate'] = 0.07
        memento.setup_memento(cell_adata, q_column='capture_rate')

        memento.create_groups(cell_adata, label_columns=['stim', 'ind'])

        memento.compute_1d_moments(cell_adata,
                                   min_perc_group=.9) # percentage of groups that satisfy the condition for a gene to be considered. 

        # Perform 1D hypothesis testing
        sample_meta = memento.get_groups(cell_adata)
        sample_meta['ind'] = sample_meta['ind'].astype('category') # make sure to not confuse ourselves in case replicate labels are numbers.
        sample_meta.head(3)

        treatment_df = sample_meta[['stim']]
        treatment_df.head(5)

        cov_df = pd.get_dummies(sample_meta['ind'].astype('category'))
        cov_df.head(3)

        memento.ht_1d_moments(
             cell_adata, 
             treatment=treatment_df,
             covariate=cov_df,
             num_boot=5000, 
             verbose=1,
             num_cpus=40)

        result_1d = memento.get_1d_ht_result(cell_adata)

        plt.scatter(result_1d.de_coef, result_1d.dv_coef, s=1)
        plt.title(f'Differential Expression for {cell_type}')
        plt.xlabel('Differential Expression Coefficient')
        plt.ylabel('Differential Variability Coefficient')
        plt.savefig(f'differential_expression_replicate_{cell_type}.pdf')
        plt.close()

        # Print the top differential mean genes (DM) and top differential variability genes (DV)
        print(f"Top DM genes for {cell_type}:")
        print(result_1d.query('de_coef > 0').sort_values('de_pval').head(10))
        
        print(f"Top DV genes for {cell_type}:")
        print(result_1d.query('dv_coef > 0 & de_coef > 0').sort_values('dv_pval').head(10))

        # Save result_1d to a txt file
        result_1d.to_csv(f'result_1d_replicate_{cell_type}.txt', sep='\t', index=False)
        
        top10_genes = result_1d.sort_values(by='dv_pval').head(10)
        print(top10_genes)
        top10_genes = top10_genes['gene']


        # Loop through research genes for 2D hypothesis testing
        for gene in top10_genes:
            gene_pairs = list(itertools.product([gene], cell_adata.var.index.tolist()))
            memento.compute_2d_moments(cell_adata, gene_pairs)
            memento.ht_2d_moments(
                cell_adata, 
                treatment=treatment_df,
                covariate=cov_df,
                num_boot=5000, 
                verbose=1,
                num_cpus=40)
            
            result_2d = memento.get_2d_ht_result(cell_adata)
            result_2d
            result_2d.sort_values('corr_pval').head(10)

            print(f"Top gene pairs for {cell_type} with {gene}:")
            print(result_2d.sort_values('corr_pval').head(5))
            result_2d.to_csv(f'result_2d_replicate_{cell_type}_{gene}.txt', sep='\t', index=False)

if __name__ == '__main__':
    main()