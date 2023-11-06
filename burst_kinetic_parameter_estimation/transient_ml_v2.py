import argparse, sys, os
import numpy as np
import ctypes
from collections import Counter
import pandas as pd
from joblib import delayed,Parallel
from scipy.optimize import minimize
def make_c_function(prec):
    c_lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'hyp1f1_arb_vector.so'))
    c_lib.calculate_probs.argtypes =  [ctypes.POINTER(ctypes.c_double), ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_ulong, ctypes.c_ulong]
    c_lib.calculate_probs.restype = None
    prec_ulong = ctypes.c_ulong(prec)
    def prepared_function(kon_double, koff_double, ksyn_double, t_double, nmax):

        kon_c_double = ctypes.c_double(kon_double)
        koff_c_double = ctypes.c_double(koff_double)
        ksyn_c_double = ctypes.c_double(ksyn_double)
        t_c_double = ctypes.c_double(t_double)
        nmax_ulong = ctypes.c_ulong(nmax)

        res = (ctypes.c_double*(nmax+1))()

        c_lib.calculate_probs(res, kon_c_double, koff_c_double, ksyn_c_double, t_c_double, nmax_ulong, prec_ulong)
        
        res_list = list(res)
        
        return res_list
    return prepared_function

def make_loglikelihood_fun(t, vals, prob_fun, set_to_zero):
    if set_to_zero:
        vals = [0 if v in set_to_zero else v for v  in vals]
    vals_dict = Counter(vals)
    nmax = max(vals)
    def loglikelihood(x):
        kon = x[0]
        koff = x[1]
        ksyn = x[2]
        
        prob_list = prob_fun(kon,koff,ksyn,t,nmax)
        for tx in set_to_zero:
            if tx < 1: continue
            try:
                prob_list[0] += prob_list[tx]
                prob_list[tx] = 1e-300
            except IndexError:
                continue
        ll_dict = dict(enumerate(np.log(np.clip(prob_list, 1e-16,1))))
        
        return -np.sum([n_obs*ll_dict[count] for count,n_obs in vals_dict.items()])
    return loglikelihood

def make_approximate_guess(bnds, grid_size, fun):
    kon_array = np.logspace(np.log10(bnds[0][0]), np.log10(bnds[0][1]), grid_size)
    koff_array = np.logspace(np.log10(bnds[1][0]), np.log10(bnds[1][1]), grid_size)
    ksyn_array = np.logspace(np.log10(bnds[2][0]), np.log10(bnds[2][1]), grid_size)
    ll_min = 1e16
    idx_min = [-1,-1,-1]
    for i, kon in enumerate(kon_array):
        for j, koff in enumerate(koff_array):
            for k, ksyn in enumerate(ksyn_array):
                if ksyn/koff < 0.1:
                    continue
                x = [kon,koff,ksyn]
                ll = fun(x)
                if ll < ll_min:
                    ll_min = ll
                    idx_min = [i,j,k]
    return [kon_array[idx_min[0]], koff_array[idx_min[1]], ksyn_array[idx_min[2]]]
def res_format(res):
    res_d = {}
    res_d['kon'] = res.x[0]
    res_d['koff'] = res.x[1]
    res_d['ksyn'] = res.x[2]
    res_d['burst_size'] = res_d['ksyn']/res_d['koff']
    res_d['ll'] = res.fun
    res_d['message'] = res.message
    res_d['success'] = res.success
    return res_d
def MaximumLikelihood(vals, t, bnds, grid_size, prob_fun, gene, x0, set_to_zero):
    if isinstance(prob_fun, int): prob_fun = make_c_function(prob_fun)
    fun = make_loglikelihood_fun(t, vals, prob_fun, set_to_zero)
    if x0 is None:
        x0 = make_approximate_guess(bnds, grid_size, fun)
    
    bnds = tuple((k/100, k*100) for k in x0)
    res = minimize(fun, x0, method='L-BFGS-B', bounds=bnds)
    res_d = res_format(res)
    return gene, res_d

def guesses(guess_df, gene, default, guess_divisor):
    if guess_df is None or gene not in guess_df.index:
        return default
    return list(guess_df.loc[gene, ['kon', 'koff', 'ksyn']]/guess_divisor)

def main():
    parser = argparse.ArgumentParser(description='Maximum likelihood inference of transient transcriptional bursting with known labelling time')
    parser.add_argument('-i', '--input', help='Input .csv file', required=True)
    parser.add_argument('-o', '--output', help='Output .csv file', required=True)
    parser.add_argument('--time', default=2, help='Labelling time of experiment (hours)', type=float)
    parser.add_argument('--threads', help='Number of threads', type=int, required=True)
    parser.add_argument('--grid_size', default=10, help='Grid size of initial approximation before L-BFGS-B algorithm')
    parser.add_argument('--prec', type=int, default=1000)
    parser.add_argument('--guesses')
    parser.add_argument('-d', '--degradation_rate', default=1)
    parser.add_argument('--set_to_zero', nargs='+', type=int, default=[])
    args = parser.parse_args()
    prob_fun = args.prec
    
    

    print('Loading file: {}'.format(args.input))
    df = pd.read_csv(args.input, index_col=0).astype(int)

    guess_df = None if args.guesses is None else pd.read_csv(args.guesses, index_col=0, sep='\t')[['kon', 'koff', 'ksyn']].dropna()
    
    try:
        degradation_rate = pd.Series(float(args.degradation_rate), index=df.index)
    except ValueError:
        degradation_rate = pd.read_csv(args.degradation_rate, index_col=0)['degradationrate']
    bnds = tuple((l/degradation_rate.mean(), u/degradation_rate.mean()) for l,u in ((1e-3,1e3),(1e-3,1e3), (1e-1, 1e3)))
    
    res = Parallel(n_jobs=args.threads, verbose=3)(delayed(MaximumLikelihood)(np.around(vals[pd.notnull(vals)]), args.time*degradation_rate[gene], bnds, args.grid_size, prob_fun, gene, guesses(guess_df, gene, None, degradation_rate[gene]), args.set_to_zero) for gene, vals in df.iterrows())

    res_dict = {r[0]: r[1] for r in res}

    res_df = pd.DataFrame.from_dict(res_dict, orient='index')
    res_df['kon'] *= degradation_rate
    res_df['koff'] *= degradation_rate
    res_df['ksyn'] *= degradation_rate

    res_df.to_csv(args.output)


if __name__ == '__main__':
    main()