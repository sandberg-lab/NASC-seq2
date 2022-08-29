import pandas, argparse, json, jsonpickle, gzip
from scipy import interpolate

def hash_dataframe(params_to_zeros_table):
    return int(pandas.util.hash_pandas_object(params_to_zeros_table).sum())

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('table_in')
    parser.add_argument('table_out', nargs='?')
    parser.add_argument('--min_expr', default=0.01, type=float)
    parser.add_argument('--max_expr', default=350.0, type=float)
    parser.add_argument('--kon', nargs=2, type=float)
    parser.add_argument('--koff', nargs=2, type=float)
    parser.add_argument('--ksyn', nargs=2, type=float)
    parser.add_argument('--json', nargs='?', const=True, help=argparse.SUPPRESS)
    o = parser.parse_args()
    
    if o.table_out is None:
        if '.trimmed.' is o.table_in: raise Exception
        
        parts = o.table_in.rsplit('.', 1)
        o.table_out = '.trimmed.'.join(parts)
    
    if o.json is True:
        o.json = o.table_out + '.json.gz'
    
    open_chance_name, close_chance_name, transcribe_chance_name, freq_adj0_name, size_adj0_name, size_cov_name = 'open_chance', 'close_chance', 'transcribe_chance', 'freq_adj0', 'size_adj0', 'size_cov'
    lookup = pandas.read_csv(o.table_in, sep='\t')
    lookup_h=lookup
    lookup_s = lookup_h[(lookup_h[freq_adj0_name]*lookup_h[size_adj0_name]<o.max_expr) & (lookup_h[freq_adj0_name]*lookup_h[size_adj0_name]>o.min_expr)]
    
    for param, limvalues in zip([open_chance_name, close_chance_name, transcribe_chance_name], [o.kon, o.koff, o.ksyn]):
        if limvalues is not None:
            lookup_s = lookup_s[(lookup_s[param]<limvalues[1]) & (lookup_s[param]>=limvalues[0])]
    
    
    params_to_zeros_table = lookup_s
    param_lookup = list(zip(params_to_zeros_table[open_chance_name], params_to_zeros_table[close_chance_name], params_to_zeros_table[transcribe_chance_name], params_to_zeros_table[freq_adj0_name], params_to_zeros_table[size_adj0_name], params_to_zeros_table[size_cov_name]))
    interpolation_points = interpolate.LinearNDInterpolator([param_lookup_entry[3:] for param_lookup_entry in param_lookup], [param_lookup_entry[:3] for param_lookup_entry in param_lookup])
    
    params_to_zeros_table.to_csv(o.table_out, sep='\t', index=False)
    
    if o.json is not None:
        with (gzip.open if o.json.endswith('.gz') else open)(o.json, 'wt') as jsonh:
            frozen = jsonpickle.encode(interpolation_points)
            json.dump({'tablehash':hash_dataframe(params_to_zeros_table), 'interpolation_points':frozen}, jsonh)