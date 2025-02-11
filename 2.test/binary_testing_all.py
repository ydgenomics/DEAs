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
from adjustText import adjust_text

@click.command()
@click.option('--input_h5ad', default='/data/work/input/memento/Pog_unsoupx_memento.h5ad', help='Path to the input h5ad file.')
@click.option('--ctrl_name', default='V2.5R2404290045', help='Control sample name.')
@click.option('--cell_type_list', default='0', help='Comma-separated list of cell types.')
@click.option('--capture_rate', default=0.07, help='Capture rate.')
@click.option('--coef_threshold', default=1.0, help='Coefficient threshold.')
@click.option('--pval_threshold', default=0.01, help='P-value threshold.')
@click.option('--top_number', default=10, help='Number of top genes to display.')
@click.option('--perform_2d_test', default='no', help='Perform 2D hypothesis testing?')

def main(input_h5ad, ctrl_name, cell_type_list, capture_rate, coef_threshold, pval_threshold, top_number, perform_2d_test):
    cell_type_list = cell_type_list.split(',')
    capture_rate = float(capture_rate)
    coef_threshold = float(coef_threshold)
    pval_threshold = float(pval_threshold)
    top_number = int(top_number)

    adata = sc.read(input_h5ad)
    print(adata)

    for cell_type in cell_type_list:
        process_cell_type(adata, cell_type, ctrl_name, capture_rate, coef_threshold, pval_threshold, top_number, perform_2d_test)

def process_cell_type(adata, cell_type, ctrl_name, capture_rate, coef_threshold, pval_threshold, top_number, perform_2d_test):
    print(f"Processing cell type: {cell_type}")
    cell_type_dir = create_directories(cell_type)
    cell_adata = filter_data(adata, cell_type, ctrl_name, capture_rate)
    result_1d = perform_1d_hypothesis_testing(cell_adata)
    plot_results(result_1d, cell_type, coef_threshold, pval_threshold)
    save_results(result_1d, cell_type)
    top_genes = get_top_genes(result_1d, top_number)
    plot_top_genes(cell_adata, top_genes, cell_type_dir, perform_2d_test, cell_type)

def create_directories(cell_type):
    cell_type_dir = os.path.join(os.getcwd(), cell_type)
    os.makedirs(cell_type_dir, exist_ok=True)
    os.chdir(cell_type_dir)
    os.makedirs('topgenes_2d_test', exist_ok=True)
    os.makedirs('topgenes_boxplot', exist_ok=True)
    return cell_type_dir

def filter_data(adata, cell_type, ctrl_name, capture_rate):
    cell_adata = adata[adata.obs['cell'] == cell_type].copy()
    cell_adata.obs['stim'] = cell_adata.obs['stim'].apply(lambda x: 0 if x == ctrl_name else 1)
    cell_adata.obs['capture_rate'] = capture_rate
    memento.setup_memento(cell_adata, q_column='capture_rate')
    memento.create_groups(cell_adata, label_columns=['stim', 'ind'])
    memento.compute_1d_moments(cell_adata, min_perc_group=.9)
    return cell_adata

def perform_1d_hypothesis_testing(cell_adata):
    sample_meta = memento.get_groups(cell_adata)
    sample_meta['ind'] = sample_meta['ind'].astype('category')
    treatment_df = sample_meta[['stim']]
    cov_df = pd.get_dummies(sample_meta['ind'].astype('category'))
    memento.ht_1d_moments(cell_adata, treatment=treatment_df, covariate=cov_df, num_boot=5000, verbose=1, num_cpus=40)
    return memento.get_1d_ht_result(cell_adata)

def plot_results(result_1d, cell_type, coef_threshold, pval_threshold):
    plt.scatter(result_1d.de_coef, result_1d.dv_coef, s=1)
    plt.title(f'Differential Expression for {cell_type}')
    plt.xlabel('Differential Expression Coefficient')
    plt.ylabel('Differential Variability Coefficient')
    plt.savefig(f'differential_expression_replicate_{cell_type}.pdf')
    plt.close()
    plot_volcano(result_1d, 'de_coef', 'de_pval', coef_threshold, pval_threshold, 'Volcano Plot for Differential Expression', f'volcano_DE_{cell_type}.pdf')
    plot_volcano(result_1d, 'dv_coef', 'dv_pval', coef_threshold, pval_threshold, 'Volcano Plot for Differential Variability', f'volcano_DV_{cell_type}.pdf')

def plot_volcano(data, coef_col, pval_col, coef_threshold, pval_threshold, title, filename):
    """
    Plots a volcano plot for differential expression analysis.

    Parameters:
    data (pd.DataFrame): DataFrame containing the differential expression data.
    coef_col (str): Column name for the coefficient (log fold change).
    pval_col (str): Column name for the p-value.
    coef_threshold (float): Threshold for the coefficient to determine significance.
    pval_threshold (float): Threshold for the p-value to determine significance.
    title (str): Title of the plot.
    filename (str): Filename to save the plot.

    Returns:
    None: The function saves the plot to the specified filename.

    The function performs the following steps:
    1. Adds a 'Significance' column to the data based on the provided thresholds.
    2. Handles zero p-values by replacing them with the minimum non-zero p-value.
    3. Calculates the negative log10 of the p-values.
    4. Creates a scatter plot with different colors for significant and non-significant points.
    5. Adds horizontal and vertical lines to indicate the thresholds.
    6. Annotates the top 10 significant genes and genes with zero p-values.
    7. Saves the plot to the specified filename.
    """
    data['Significance'] = np.select(
        [(data[coef_col] > coef_threshold) & (data[pval_col] < pval_threshold),
         (data[coef_col] < -coef_threshold) & (data[pval_col] < pval_threshold),
         (data[pval_col] >= pval_threshold)],
        ['Significant Up', 'Significant Down', 'Not Significant'],
        default='Not Significant'
    )
    zero_pval_col_genes = data.loc[data[pval_col] == 0]
    if not zero_pval_col_genes.empty:
        min_nonzero_pval = data.loc[data[pval_col] > 0, pval_col].min()
        if pd.isna(min_nonzero_pval):
            min_nonzero_pval = 1e-10
        data.loc[data[pval_col] == 0, pval_col] = min_nonzero_pval
    data['neg_log_pval'] = -np.log10(data[pval_col].replace(0, np.nan).fillna(data[pval_col].min()))
    zero_pval_col_genes = data[data['gene'].isin(zero_pval_col_genes['gene'])]
    print("zero pval genes:")
    print(zero_pval_col_genes)
    plt.figure(figsize=(10, 6), dpi=1200)
    sns.scatterplot(x=coef_col, y='neg_log_pval', hue='Significance', data=data,
                    palette={'Significant Up': 'red', 'Significant Down': 'blue', 'Not Significant': 'gray'})
    plt.axhline(y=-np.log10(pval_threshold), color='black', linestyle='--', label=f'p={pval_threshold}')
    plt.axvline(x=coef_threshold, color='black', linestyle='--', label=f'{coef_col}={coef_threshold}', alpha=0.5)
    plt.axvline(x=-coef_threshold, color='black', linestyle='--', alpha=0.5)
    significant_genes = data[(data['Significance'] == 'Significant Up') | (data['Significance'] == 'Significant Down')]
    top10_significant_genes = significant_genes.nsmallest(10, pval_col)
    texts2 = []
    print(zero_pval_col_genes)
    if not zero_pval_col_genes.empty:
        for _, row in zero_pval_col_genes.iterrows():
            texts2.append(plt.text(row[coef_col], row['neg_log_pval'], row['gene'], fontsize=8, ha='center', va='center', color='green'))
        adjust_text(texts2, arrowprops=dict(arrowstyle='->', color='gray'))
    zero_pval_genes_set = set(zero_pval_col_genes['gene'])
    top10_significant_genes_set = set(top10_significant_genes['gene'])
    non_intersect_genes = top10_significant_genes_set - zero_pval_genes_set
    non_intersect_top10_significant_genes = top10_significant_genes[top10_significant_genes['gene'].isin(non_intersect_genes)]
    texts1 = []
    for _, row in non_intersect_top10_significant_genes.iterrows():
        texts1.append(plt.text(row[coef_col], row['neg_log_pval'], row['gene'], fontsize=8, ha='center', va='center', color='black'))
    adjust_text(texts1, arrowprops=dict(arrowstyle='->', color='gray'))
    plt.xlabel(coef_col)
    plt.ylabel('-Log10(P-value)')
    plt.title(title)
    plt.legend(title='Significance')
    plt.savefig(filename)
    plt.close()

def save_results(result_1d, cell_type):
    result_1d.to_csv(f'result_1d_replicate_{cell_type}.txt', sep='\t', index=False)

def get_top_genes(result_1d, top_number):
    return result_1d.sort_values(by='dv_pval').head(top_number)['gene']

def plot_top_genes(cell_adata, top_genes, cell_type_dir, perform_2d_test, cell_type):
    mean, var, counts = memento.get_1d_moments(cell_adata)
    for gene in top_genes:
        plot_gene_moments(mean, var, gene, cell_type_dir)
        if perform_2d_test == 'yes':
            perform_2d_hypothesis_testing(cell_adata, gene, cell_type_dir, cell_type)
    os.chdir(cell_type_dir)
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

def perform_2d_hypothesis_testing(cell_adata, gene, cell_type_dir, cell_type):
    os.chdir(os.path.join(cell_type_dir, 'topgenes_2d_test'))
    sample_meta = memento.get_groups(cell_adata)
    sample_meta['ind'] = sample_meta['ind'].astype('category')
    treatment_df = sample_meta[['stim']]
    cov_df = pd.get_dummies(sample_meta['ind'].astype('category'))
    #result_2d.to_csv(f'result_2d_replicate_{cell_adata.obs["cell"].unique()[0]}_{gene}.txt', sep='\t', index=False)
    gene_pairs = list(itertools.product([gene], cell_adata.var.index.tolist()))
    memento.compute_2d_moments(cell_adata, gene_pairs)
    memento.ht_2d_moments(cell_adata, treatment=treatment_df, covariate=cov_df, num_boot=5000, verbose=1, num_cpus=40)
    result_2d = memento.get_2d_ht_result(cell_adata)
    result_2d.to_csv(f'result_2d_replicate_{cell_type}_{gene}.txt', sep='\t', index=False)

if __name__ == '__main__':
    main()
