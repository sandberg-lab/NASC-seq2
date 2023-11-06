import collections, itertools
import numpy as np


def make_c_function(prec):
    import ctypes
    c_lib = ctypes.CDLL('/home/danielr/cdk9as/src/hyp1f1_arb-main/hyp1f1_arb_vector.so')
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

        return list(res)
    return prepared_function

Burstnumbers = collections.namedtuple('Burstnumbers', ['freq_adj0', 'size_adj0', 'size_cov', 'open_chance', 'close_chance', 'transcribe_chance', 'degrade_chance', 'loss_chance', 'burst_size', 'mean_occupancy', 'switching_rate'])

def calc_numbers(sim_output_list, open_chance, close_chance, transcribe_chance, degrade_chance):
    fraction_nonzero = 1 - sim_output_list[0]
    size_raw = sum(tx*prob for tx, prob in enumerate(sim_output_list))/fraction_nonzero
    size_cov = (sum((tx - size_raw)**2 * prob for tx, prob in enumerate(sim_output_list) if tx > 0)/fraction_nonzero)**0.5/size_raw # coefficient of variance, excluding zeros
    return Burstnumbers(fraction_nonzero, size_raw, size_cov, open_chance, close_chance, transcribe_chance, degrade_chance, 0, transcribe_chance/close_chance, open_chance/(open_chance+close_chance), 1/(open_chance+close_chance))


def sim_and_calc_nums(open_chance, close_chance, transcribe_chance, degrade_chance, time, maxexpr, prob_fun):
    try:
        params = open_chance/degrade_chance, close_chance/degrade_chance, transcribe_chance/degrade_chance, time*degrade_chance, maxexpr
        probs = prob_fun(*params)
        for tx, prob in enumerate(probs):
            if prob < 0 or prob > 1:
                if tx-10 < maxexpr/10:
                    return None
                else:
                    probs = list(probs)[:tx-10]
        return calc_numbers(probs, open_chance, close_chance, transcribe_chance, degrade_chance)._asdict()
    except OverflowError:
        return None

def sim_and_calc_param_list(param_list, prec, jobnum):
    prob_fun = make_c_function(prec)
    ret = [sim_and_calc_nums(*params, prob_fun) for params in param_list]
    return [r for r in ret if r is not None]

def intertwine(*args):
    return (item for pair in itertools.zip_longest(*args) for item in pair if item is not None)

def chances(option):
    minimum, maximum, count = option
    if count == 1:
        return (minimum, )
    elif count == 2:
        return (minimum, maximum)
    else:
        step = (np.log2(maximum)-np.log2(minimum))/(count-1)
        ordered = list(np.exp2(np.arange(np.log2(minimum), np.log2(maximum)+step/2, step)))
        return intertwine(ordered[:len(ordered)//4], ordered[len(ordered)//4:len(ordered)//2], reversed(ordered[len(ordered)//2:3*len(ordered)//4]), reversed(ordered[3*len(ordered)//4:]))

if '__main__' == __name__:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--plot_out')
    parser.add_argument('--open_rate', type=float, nargs=3, help='kon', required=True, metavar=('min', 'max', 'count'))
    parser.add_argument('--close_rate', type=float, nargs=3, help='koff', default=[3, 3, 1], metavar=('min', 'max', 'count'))
    parser.add_argument('--transcribe_rate', type=float, nargs=3, help='ksyn', required=True, metavar=('min', 'max', 'count'))
    parser.add_argument('--degrade_rate', type=float, nargs=3, default=[0.015, 0.015, 1], metavar=('min', 'max', 'count'))
    parser.add_argument('--time', type=float, nargs=3, default=[1, 1, 1], metavar=('min', 'max', 'count'))
    parser.add_argument('--table_out', required=True)
    parser.add_argument('--proc', type=int, default=1)
    parser.add_argument('--maxexpr', type=int, default=1000)
    parser.add_argument('--prec', type=int, default=10000)
    o = parser.parse_args()
    
    import tqdm, pandas, math
    from concurrent import futures
    
    
    param_list = []
    for open_chance in chances(o.open_rate):
        for close_chance in chances(o.close_rate):
            for transcribe_chance in chances(o.transcribe_rate):
                for degrade_chance in chances(o.degrade_rate):
                    for time in chances(o.time):
                        param_list.append((open_chance, close_chance, transcribe_chance, degrade_chance, time, o.maxexpr))
    
    if o.proc == 1:
        prob_fun = make_c_function(o.prec)
        before_nums_list = [sim_and_calc_nums(*params, prob_fun) for params in tqdm.tqdm(param_list, colour='#cccccc')]
    else:
        perproc = int(math.ceil(len(param_list)/2/o.proc))
        
        pool = futures.ProcessPoolExecutor(o.proc)
        jobs = []
        while param_list:
            list_len_taken = min(len(param_list), perproc)
            jobs.append(pool.submit(sim_and_calc_param_list, param_list[:list_len_taken], o.prec, len(jobs)))
            param_list = param_list[list_len_taken:]
        before_nums_list = []
        for job in tqdm.tqdm(jobs, position=0):
            before_nums_list.extend(job.result())
    
    table_out = pandas.DataFrame()
    for nums_list in (before_nums_list, ):
        for nums in nums_list:
            table_out = table_out.append(nums, ignore_index=True)
    table_out.to_csv(o.table_out, sep='\t', index=False)
    
    if o.plot_out is not None:
        from matplotlib import pyplot
        import seaborn
        pyplot.clf()
        seaborn.scatterplot(x='freq_adj0', y='size_adj0', hue='size_cov', data=table_out)
        pyplot.xlim(0, 1)
        pyplot.ylim(bottom=1)
        pyplot.savefig(o.plot_out)