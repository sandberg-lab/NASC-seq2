import pandas, numpy, os
from scipy import interpolate

def loadlist(path):
    with open(path, 'r') as infh:
        return [l.strip() for l in infh]

def binary(value):
    return int(value > 0)

def do_some_permutations(func, c, lena, n):
    import random
    arr = []
    while len(arr) < n:
        random.shuffle(c)
        arr.append(abs(func(c[:lena]) - func(c[lena:])))
    return arr

def bootstrap_ci(func, a, n_ctrl):
    import random
    arr =[]
    while len(arr) <= n_ctrl:
        try:
            #c = random.choices(a, k=len(a))
            c = a.sample(frac=1, replace=True)
        except AttributeError:
            raise
        arr.append(func(c))
    noNA = [v for v in arr if str(v) != 'nan']
    
    if len(noNA) < 10:
        return float('nan'), float('nan')
    
    conf_width_values = numpy.percentile(noNA, [2.5, 97.5])
    
    return min(conf_width_values), max(conf_width_values)

def kdtree_query_and_interpolation(kdtree, point, index, param_lookup, interpolation_points, loginterpolate):
    if kdtree is None: raise ValueError
    import warnings
    warnings.filterwarnings("ignore")
    try:
        v = interpolation_points(point)[index]
        if 'nan' in str(v): raise Exception 
        ret = v
    except:
        d, closest_lookup_index = kdtree.query(point)
        ret = param_lookup[closest_lookup_index][index]
    if loginterpolate:
        return 2**ret
    else:
        return ret

def calc_size_cov(value_series, size_adj0, add_cov):
    cov = numpy.std([v for v in value_series if v > 0])/size_adj0
    return max(0, cov+add_cov)

def burstnum_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series):
    freq_adj0 = numpy.mean([binary(v) for v in value_series])
    if freq_adj0 == 0: return 0
    if freq_adj0 == 1: return float('inf')
    size_adj0 = value_series.mean() / freq_adj0
    size_cov = calc_size_cov(value_series, size_adj0, add_cov)
    return kdtree_query_and_interpolation(freq_size_adj0_kdtree, (freq_adj0, size_adj0, size_cov), 0, param_lookup, interpolation_points, loginterpolate)

def size_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series):
    return ksyn_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series)/closing_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series)

def ksyn_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series):
    freq_adj0 = numpy.mean([binary(v) for v in value_series])
    if freq_adj0 == 0: return float('nan')
    if freq_adj0 == 1: return float('nan')
    size_adj0 = value_series.mean() / freq_adj0
    size_cov = calc_size_cov(value_series, size_adj0, add_cov)
    return kdtree_query_and_interpolation(freq_size_adj0_kdtree, (freq_adj0, size_adj0, size_cov), 1, param_lookup, interpolation_points, loginterpolate)

def closing_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series):
    freq_adj0 = numpy.mean([binary(v) for v in value_series])
    if freq_adj0 == 0: return float('nan')
    if freq_adj0 == 1: return float('nan')
    size_adj0 = value_series.mean() / freq_adj0
    size_cov = calc_size_cov(value_series, size_adj0, add_cov)
    return kdtree_query_and_interpolation(freq_size_adj0_kdtree, (freq_adj0, size_adj0, size_cov), 2, param_lookup, interpolation_points, loginterpolate)

def calc_nums(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series):
    freq_adj0 = value_series.apply(binary).mean()
    if freq_adj0 == 0:
        return 0, 0, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')
    size_adj0 = value_series.mean() / freq_adj0
    size_cov = calc_size_cov(value_series, size_adj0, add_cov)
    if freq_adj0 == 1:
        return len(value_series), 1, size_adj0, size_cov, float('inf'), float('nan'), float('nan'), float('nan')
    burstnum = burstnum_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series)
    closing =   closing_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series)
    ksyn =         ksyn_from_nonzeros(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, value_series)
    burstsize =    ksyn/closing
    return value_series.apply(binary).sum(), freq_adj0, size_adj0, size_cov, burstnum, closing, ksyn, burstsize

def one_gene(gene, df1_gene, freq_size_adj0_kdtree, param_lookup, add_cov, permutations, interpolation_points, loginterpolate, proc):
    import warnings, functools
    #warnings.filterwarnings("ignore")
    
    if freq_size_adj0_kdtree is None:
        freq_size_adj0_kdtree, param_lookup, interpolation_points = build_KD_tree(param_lookup, loginterpolate)
    
    
    row = pandas.Series(name=gene)
    row['binary_sum'], row["binary_avg"], row["nonzero_avg"], row["nonzero_cov"], row[burstnum_name], row[koff_name], row[ksyn_name], row[burstsize_name] = calc_nums(freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov, df1_gene)
    
    for param, lookup_func in [(burstnum_name, burstnum_from_nonzeros), (koff_name, closing_from_nonzeros), (ksyn_name, ksyn_from_nonzeros), (burstsize_name, size_from_nonzeros)]:
        row[param+'_95%conf_low'], row[param+'_95%conf_high'] = bootstrap_ci(functools.partial(lookup_func, freq_size_adj0_kdtree, param_lookup, interpolation_points, loginterpolate, add_cov), df1_gene, permutations)
    
    return row


def build_KD_tree(params_to_zeros, loginterpolate):
    params_to_zeros_table = pandas.read_table(params_to_zeros).dropna().drop_duplicates(subset=[freq_adj0_name, size_adj0_name, size_cov_name], keep=False)
    try:
        freq_size_adj0_kdtree = spatial.KDTree(params_to_zeros_table[[freq_adj0_name, size_adj0_name, size_cov_name]])
    except RecursionError:
        warn('Could not create a KDtree due to a RecursionError')
        freq_size_adj0_kdtree = 1
    if loginterpolate:
        param_lookup = list(zip(numpy.log2(params_to_zeros_table[open_chance_name]), numpy.log2(params_to_zeros_table[transcribe_chance_name]), numpy.log2(params_to_zeros_table[close_chance_name]), numpy.log2(params_to_zeros_table[burst_size_name]), params_to_zeros_table[freq_adj0_name], params_to_zeros_table[size_adj0_name], params_to_zeros_table[size_cov_name]))
    else:
        param_lookup = list(zip(params_to_zeros_table[open_chance_name], params_to_zeros_table[transcribe_chance_name], params_to_zeros_table[close_chance_name], params_to_zeros_table[burst_size_name], params_to_zeros_table[freq_adj0_name], params_to_zeros_table[size_adj0_name], params_to_zeros_table[size_cov_name]))
    interpolation_points = interpolate.LinearNDInterpolator([param_lookup_entry[4:] for param_lookup_entry in param_lookup], [param_lookup_entry[:4] for param_lookup_entry in param_lookup])
    return freq_size_adj0_kdtree, param_lookup, interpolation_points

if '__main__' == __name__:
    import tqdm
    import collections, argparse
    from scipy import spatial
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--params_to_zeros", required=True)
    parser.add_argument('--max_genes', type=int)
    parser.add_argument('-bootstraps', '--permutations', type=int, default=1000, dest='permutations')
    parser.add_argument('--proc', nargs=2, type=int, default=[10, 1])
    parser.add_argument("csv1", metavar='csv_G1')
    parser.add_argument("table_out")
    parser.add_argument('--reuse_lookups', action='store_true')
    parser.add_argument('--parameter_names', default=['freq_adj0', 'size_adj0', 'size_cov', 'open_chance', 'close_chance', 'transcribe_chance', 'burst_size', 'kon', 'koff', 'ksyn', 'burstsize'], metavar=('input_nonzero_fraction', 'input_nonzero_average', 'input_nonzero_cov', 'input1', 'input2', 'input3', 'input4', 'output1', 'output2', 'output3', 'output4'), nargs=11)
    parser.add_argument('--simcells', action='store_true')
    parser.add_argument('--loginterpolate', action='store_true')
    parser.add_argument('--add_cov', type=float, default=0)
    o = parser.parse_args()
    load_lookups_per_gene = not o.reuse_lookups
    
    if os.path.exists(o.table_out): raise Exception
    
    global freq_adj0_name, size_adj0_name, size_cov_name, open_chance_name, close_chance_name, transcribe_chance_name, burst_size_name, burstnum_name, koff_name, ksyn_name, burstsize_name
    freq_adj0_name, size_adj0_name, size_cov_name, open_chance_name, close_chance_name, transcribe_chance_name, burst_size_name, burstnum_name, koff_name, ksyn_name, burstsize_name = o.parameter_names
    
    # load data
    if o.simcells:
        df1 = pandas.read_table(o.csv1, nrows=o.max_genes)
        df1.index = ['simgene%i_kon%.3f_koff%.3f_ksyn%.3f_expr%.3f'%(i, kon, koff, ksyn, expr) for i, (kon, koff, ksyn, expr) in enumerate(zip(df1['open_chance'], df1['close_chance'], df1['transcribe_chance'], df1['freq_adj0'] * df1['size_adj0']))]
        df1 = df1[[col for col in df1.columns if col.startswith('simcell')]]
        df1 = df1.T
    else:
        df1 = pandas.read_csv(o.csv1, index_col=0, nrows=o.max_genes).T
    
    
    if load_lookups_per_gene:
        freq_size_adj0_kdtree = None
        param_lookup = o.params_to_zeros
        interpolation_points = None
    else:
        freq_size_adj0_kdtree, param_lookup, interpolation_points = build_KD_tree(o.params_to_zeros, o.loginterpolate)
    
    # test on frequency of non-zero
    from concurrent import futures
    pool = futures.ProcessPoolExecutor(o.proc[0])
    jobs = []
    for gene in df1:
        jobs.append(pool.submit(one_gene, gene, df1[gene], freq_size_adj0_kdtree, param_lookup, o.add_cov, o.permutations, interpolation_points, o.loginterpolate, o.proc[1]))
    table_out = pandas.DataFrame([job.result() for job in tqdm.tqdm(jobs)])
    
    table_out.to_csv(o.table_out, sep='\t')
    
    