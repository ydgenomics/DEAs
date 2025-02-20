import scanpy as sc
import click

@click.command()
@click.option('--input_h5ad', default=None, help='Path to h5ad file')
@click.option('--pre_ind', default='pre_ind', help='Column name for pre_ind')
@click.option('--pre_stim', default='pre_stim', help='Column name for pre_stim')
@click.option('--pre_cell', default='pre_cell', help='Column name for pre_cell')
@click.option('--out_h5ad', default=None, help='Path to output h5ad file')
def process_h5ad(input_h5ad, pre_ind, pre_stim, pre_cell, out_h5ad):
    # Read h5ad file
    adata = sc.read_h5ad(input_h5ad)
    
    # Check and assign raw counts
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()

    # Check adata.X whether is csr_matrix and change into csr_matrix
    if not isinstance(adata.X, csr_matrix):
        adata.X = csr_matrix(adata.X)

    # Check and assign observation values
    if 'ind' not in adata.obs.columns and pre_ind in adata.obs.columns:
        adata.obs['ind'] = adata.obs[pre_ind]
    
    if 'stim' not in adata.obs.columns and pre_stim in adata.obs.columns:
        adata.obs['stim'] = adata.obs[pre_stim]
        
    if 'cell' not in adata.obs.columns and pre_cell in adata.obs.columns:
        adata.obs['cell'] = adata.obs[pre_cell]
    
    print(adata)

    adata.write(out_h5ad, compression='gzip')

if __name__ == "__main__":
    process_h5ad()