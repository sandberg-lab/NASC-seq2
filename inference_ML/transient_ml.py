import argparse
import numpy as np
import ctypes
from collections import Counter
import pandas as pd
from joblib import delayed,Parallel
from scipy.optimize import minimize
def make_c_function(prec):
    c_lib = ctypes.CDLL('/home/antonl/programs/hyp_new_git/hyp1f1_arb/hyp1f1_arb_vector.so')
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

def make_loglikelihood_fun(t, vals, prob_fun):
    vals_dict = Counter(vals)
    nmax = max(vals)
    def loglikelihood(x):
        kon = x[0]
        koff = x[1]
        ksyn = x[2]
        
        prob_list = prob_fun(kon,koff,ksyn,t,nmax)
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
def MaximumLikelihood(vals, t, bnds, grid_size, prob_fun, gene):
    fun = make_loglikelihood_fun(t, vals, prob_fun)
    x0 = make_approximate_guess(bnds, grid_size, fun)
    res = minimize(fun, x0, method='L-BFGS-B', bounds=bnds)
    res_d = res_format(res)
    return gene, res_d

def main():
    parser = argparse.ArgumentParser(description='Maximum likelihood inference of transient transcriptional bursting with known labelling time')
    parser.add_argument('-i', '--input', help='Input .csv file')
    parser.add_argument('-o', '--output', help='Output .csv file')
    parser.add_argument('--time', default=2, help='Labelling time of experiment (hours)')
    parser.add_argument('--threads', help='Number of threads')
    parser.add_argument('--grid_size', default=10, help='Grid size of initial approximation before L-BFGS-B algorithm')
    bnds = ((1e-3,1e3),(1e-3,1e3), (1e-1, 1e3))
    prob_fun = make_c_function(256)
    args = parser.parse_args()

    print('Loading file: {}'.format(args.input))
    df = pd.read_csv(args.input, index_col=0)

    res = Parallel(n_jobs=args.threads, verbose=3)(delayed(MaximumLikelihood(np.around(vals[pd.notnull(vals)]), args.time, bnds, args.grid_size, prob_fun, gene)) for gene, vals in df.iterrows())

    res_dict = {r[0]: r[1] for r in res}

    res_df = pd.DataFrame.from_dict(res_dict, orient='index')

    res_df.to_csv(args.output)


if __name__ == '__main__':
    main()