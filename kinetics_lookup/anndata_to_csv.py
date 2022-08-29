import anndata, argparse, sys, random, pandas

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('anndata')
    parser.add_argument('layer')
    parser.add_argument('--rand_rows', type=int)
    parser.add_argument('--rand_cols', type=int)
    parser.add_argument('--rand_seed', type=int, default=0)
    parser.add_argument('--transpose', action='store_true')
    parser.add_argument('--invert_rand_rows', const=True, nargs='?', type=int)
    parser.add_argument('--invert_rand_cols', const=True, nargs='?', type=int)
    o = parser.parse_args()
    
    random.seed(o.rand_seed)
    
    adata = anndata.read_h5ad(o.anndata, backed='r')
    try:
        sparse_matrix = adata.layers[o.layer] # type is scipy.sparse.csr.csr_matrix
        rows = list(range(len(adata.obs.index)))
        cols = list(range(len(adata.var.index)))
        if o.rand_rows is not None:
            rows = random.sample(rows, o.rand_rows)
            if o.invert_rand_rows is not None:
                rows = [row for row in list(range(len(adata.obs.index))) if row not in rows]
                if o.invert_rand_rows is not True:
                    rows = random.sample(rows, o.invert_rand_rows)
            sparse_matrix = sparse_matrix[rows, :]
        if o.rand_cols is not None:
            cols = random.sample(cols, o.rand_cols)
            if o.invert_rand_cols is not None:
                cols = [col for col in list(range(len(adata.var.index))) if col not in cols]
                if o.invert_rand_cols is not True:
                    cols = random.sample(cols, o.invert_rand_cols)
            sparse_matrix = sparse_matrix[:, cols]
        
        df = pandas.DataFrame(sparse_matrix.A, columns=adata.var.index[cols], index=adata.obs.index[rows])
        if o.transpose:
            df = df.T
        df.to_csv(sys.stdout)
    except BrokenPipeError:
        pass
    finally:
        adata.file.close()