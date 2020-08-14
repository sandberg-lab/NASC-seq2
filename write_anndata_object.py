import anndata
import h5py
import pandas as pd
import numpy as np
from scipy import sparse
import tensorflow as tf
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract molecule information and store in an hdf5 file')
    parser.add_argument('-h5','--hdf5',metavar='input', type=str, help='.h5 file to parse')
    parser.add_argument('-ann', '--anndata', metavar='output', type=int, default=1, help='anndata file to write')
    args = parser.parse_args()
    h5file = args.hdf5
    anndata_file = parser.anndata
    h5file = h5py.File(h5file, 'r')
    print('Parsing hdf5 file: {}'.format(h5file))
    all_dicts = {'total': {}, 'new': {}, 'old': {}, 'spliced': {}, 'unspliced': {}}
    for gene, gene_grp in h5file['genes'].items():
        cells = [cell.decode('utf-8') for cell in gene_grp['cell'][:]]
        new_tensor = tf.ragged.constant(gene_grp['htest'][:])
        intronic_reads = tf.ragged.constant(gene_grp['ir'][:])
    
        new = tf.math.reduce_sum(new_tensor, axis=1).numpy()
        old = tf.math.reduce_sum(1-new_tensor, axis=1).numpy()
    
        intronic_any = tf.clip_by_value(intronic_reads,0,1)
        unspliced = tf.math.reduce_sum(intronic_any,axis=1).numpy()
        spliced = tf.math.reduce_sum(1-intronic_any,axis=1).numpy()
     
        total = new+old
    
        all_dicts['total'][gene] = pd.Series(total, index=cells)
        all_dicts['new'][gene] = pd.Series(new, index=cells)
        all_dicts['old'][gene] = pd.Series(old, index=cells)
        all_dicts['unspliced'][gene] = pd.Series(unspliced, index=cells)
        all_dicts['spliced'][gene] = pd.Series(spliced, index=cells)
    all_dfs = {}
    for k, d in all_dicts.items():
        all_dfs[k] = pd.DataFrame(d).fillna(0)
    X = all_dfs['total'].values
    X = sparse.csr_matrix(X)
    ad = anndata.AnnData(X)
    ad.obs_names = all_dfs['total'].index.values
    ad.var_names = all_dfs['total'].columns.values
    for k,df in all_dfs.items():
        if k == 'total':
            continue
        ad.layers[k] = sparse.csr_matrix(df.values)

    print('Writing to anndata file: '.format(anndata_file))
    ad.write(anndata_file)
