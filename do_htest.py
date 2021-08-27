import h5py
import numpy as np
from scipy.stats import binom
import argparse
from multiprocessing import Process, Manager
from joblib import delayed, Parallel
import os
from itertools import accumulate
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'


def calc_power(n, criteria, p_c, p_e):
    r = np.arange(n+1)
    t = np.isclose(binom.pmf(r, n, p_e)/binom.pmf(r, n, p_c), criteria['c'])
    i = np.argmax(t > 0)
    pmf_array = binom.pmf(r, n, p_c)
    lt_prob = np.sum(pmf_array[i+1:])
    eq_prob = pmf_array[i]
    power = criteria['q']*eq_prob + lt_prob
    return power


def do_htest(ll, c_array, q_array):
    close_to_c = np.isclose(ll, c_array)
    passes_random = np.random.uniform(size=len(q_array)) < q_array
    pass_but_close = close_to_c*passes_random
    passes_ll = ll < c_array
    pass_not_close = passes_ll*(~close_to_c)
    return pass_not_close + pass_but_close


def get_c_and_p(n, p_c, p_e, a):
    r = np.arange(n+1)
    sf_array = binom.sf(r, n, p_e) < a
    pmf_array = binom.pmf(r, n, p_e)
    ll_ratio_array = pmf_array/binom.pmf(r, n, p_c)
    i = np.argmax(sf_array > 0)
    lt_prob = np.sum(pmf_array[i+1:])
    eq_prob = pmf_array[i]
    p = (a-lt_prob)/eq_prob
    c = ll_ratio_array[i]
    return {'c': c, 'q': p}


def cell_new_htest(cell_id, h5file, a, q):
    fd = h5py.File(h5file, 'r')
    grp = fd['cells/{}'.format(cell_id)]
    p_c = grp.attrs['p_c']
    p_e = grp.attrs['p_e']
    if p_c < p_e:
        p_c = p_e + 1e-5
    gene_array = grp['gene'][:]
    sc = grp['sc'][:]
    len_to_split = [len(c) for c in sc]
    sc_array = np.concatenate(sc)
    tc_array = np.concatenate(grp['tc'][:])
    fd.close()

    pc_log = binom.logpmf(sc_array, tc_array, p_c)
    pe_log = binom.logpmf(sc_array, tc_array, p_e)
    rel_ll = np.exp(pe_log - pc_log)
    tc_set = set(tc_array)
    criteria_dict = {n: get_c_and_p(n, p_c, p_e, a) for n in tc_set}
    criteria_list = [(criteria_dict[tc]['c'], criteria_dict[tc]['q']) for tc in tc_array]
    power_dict = {n: calc_power(n, criteria_dict[n], p_c,p_e) for n in tc_set}
    power_array = np.array([power_dict[tc] for tc in tc_array])
    c_array = np.array([m[0] for m in criteria_list])
    q_array = np.array([m[1] for m in criteria_list])
    htest_array = do_htest(rel_ll, c_array, q_array)
    res_dict = {}
    res_dict['gene'] = gene_array
    res_dict['htest'] = split_list(htest_array, len_to_split)
    res_dict['rel_ll'] = split_list(rel_ll, len_to_split)
    res_dict['c_htest'] = split_list(c_array, len_to_split)
    res_dict['q_htest'] = split_list(q_array, len_to_split)
    res_dict['power'] = split_list(power_array, len_to_split)
    fd.close()
    q.put((cell_id, res_dict))
    return cell_id


def split_list(l, len_to_split):
    return [l[x - y: x] for x, y in zip(accumulate(len_to_split), len_to_split)]


def make_write_function(h5file):
    def write_hdf5(q):
        full_gene_dict = {}
        while True:
            cell_id, cell_dict = q.get()
            if cell_dict is None:
                break
            for gene, htest, rel_ll, c_htest, q_htest, power in zip(*cell_dict.values()):
                gene = gene.decode('utf-8')
                if gene in full_gene_dict:
                    full_gene_dict[gene][cell_id] = {'htest': htest, 'rel_ll': rel_ll, 'c_htest': c_htest, 'q_htest': q_htest, 'power': power}
                else:
                    full_gene_dict[gene] = {}
                    full_gene_dict[gene][cell_id] = {'htest': htest, 'rel_ll': rel_ll, 'c_htest': c_htest, 'q_htest': q_htest, 'power': power}
            del cell_dict
            q.task_done()
        data_to_write = {}
        hdf5_file = h5py.File(h5file, 'a')
        genes_grp = hdf5_file['genes']
        for gene, gene_grp in genes_grp.items():
            data_to_write = {'htest': [], 'rel_ll': [], 'c_htest': [], 'q_htest': [], 'power': []}
            specific_gene_dict = full_gene_dict[gene]
            for cell_id in gene_grp['cell'][:]:
                cell_id = cell_id.decode('utf-8')
                for k in data_to_write.keys():
                    data_to_write[k].append(specific_gene_dict[cell_id][k])
            for k, v in data_to_write.items():
                if k == 'htest':
                    dset = gene_grp.create_dataset(k, (len(v),), dtype=h5py.vlen_dtype(np.dtype('bool')))
                else:
                    dset = gene_grp.create_dataset(k, (len(v),), dtype=h5py.vlen_dtype(np.dtype('float32')))
                for i, arr in enumerate(v):
                    dset[i] = arr
        hdf5_file.close()
        q.task_done()
    return write_hdf5


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Do hypothesis test for new RNA')
    parser.add_argument('-h5', '--hdf5', metavar='input', type=str, help='.h5 file to process')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    h5file = args.hdf5
    threads = args.threads
    print('Doing hypothesis test for new RNA')
    fd = h5py.File(h5file, 'r')
    cell_list = list(fd['cells'].keys())
    fd.close()
    m = Manager()
    q = m.JoinableQueue()
    p = Process(target=make_write_function(h5file), args=(q,))
    p.start()
    res = Parallel(n_jobs=threads, verbose=3, backend='loky')(delayed(cell_new_htest)(cell_id, h5file, 0.01, q) for cell_id in cell_list)
    q.put((None, None))
    print('Done, writing results')
    p.join()
