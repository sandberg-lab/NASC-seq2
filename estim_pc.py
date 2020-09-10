import h5py
from joblib import delayed,Parallel
import numpy as np
from multiprocessing import Process, Manager
from scipy import sparse
import argparse
import os
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
from scipy.stats import binom
def createMkn(Akn,p_e): #Left out from Akn
    M=np.zeros(Akn.shape)
    for n in range(Akn.shape[1]):
        for k in range(Akn.shape[0]):
            Ekn= np.sum(Akn[(k+1):,n])*binom.pmf(k,n,p_e)
            if Ekn > 0.01*Akn[k,n]:
                M[k,n]=1
    return M

def EstepAkn(Akn,Mkn,p_c): # Alters Akn in place - modify initial step
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            if Mkn[k,n]==1:
                num=0
                denom=0
                for kp in range(Mkn.shape[0]):
                    if Mkn[kp,n]==1:
                        num = num + binom.pmf(k,n,p_c)*Akn[kp,n]
                        denom = denom + binom.pmf(kp,n,p_c)
                Akn[k,n]=num/denom
    return Akn

def MstepP_c(Akn):
    num=0
    denom=0
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            num=num+k*Akn[k,n]
            denom=denom+n*Akn[k,n]
    p_c=num/denom
    return p_c

#Bisection-search for p_c  #Move code above into method to get local namespace
def estimateP_c(Akn,p_e,cell_id):
    l=p_e+1e-7
    r=1
    p_c0=(l+r)/2

    
    Akn0=Akn

    p_c=p_c0
    Mkn=createMkn(Akn0,p_e)
    Akn=Akn0

    while r-l >= 10e-8:
        Akn=EstepAkn(Akn,Mkn,p_c)
        p_c_old=p_c
        p_c=MstepP_c(Akn)
        if p_c < p_c_old:
            r=p_c
        else:
            l=p_c
    return p_c, p_e, cell_id

def get_cell_dicts(genes,h5file,q):
    cell_dict = {}
    ragged_h5_local = h5py.File(h5file,'r')
    for gene in genes:
        gene_grp = ragged_h5_local['genes/{}'.format(gene)]
        if gene_grp.attrs['strand'] == '+':
            sc_array = gene_grp['sc/tC'][:]
            tc_array = gene_grp['tc/t'][:]
        else:
            sc_array = gene_grp['sc/aG'][:]
            tc_array = gene_grp['tc/a'][:]
        gA_array = gene_grp['sc/gA'][:]
        cT_array = gene_grp['sc/cT'][:]
        g_array = gene_grp['tc/g'][:]
        c_array = gene_grp['tc/c'][:]
        cells = gene_grp['cell'][:]
        for cell,sc,tc,gA,cT,g,c in zip(cells,sc_array,tc_array,gA_array,cT_array,g_array,c_array):
            cell = cell.decode('utf-8')
            if cell in cell_dict:
                cell_dict[cell]['gene'].append(gene)
                cell_dict[cell]['sc'].append(sc)
                cell_dict[cell]['tc'].append(tc)
                cell_dict[cell]['gA'] += np.sum(gA)
                cell_dict[cell]['cT'] += np.sum(cT)
                cell_dict[cell]['g'] += np.sum(g)
                cell_dict[cell]['c'] += np.sum(c)
            else:
                cell_dict[cell] = {}
                cell_dict[cell]['gene'] = [gene]
                cell_dict[cell]['sc'] = [sc]
                cell_dict[cell]['tc'] = [tc]
                cell_dict[cell]['gA'] = np.sum(gA)
                cell_dict[cell]['cT'] = np.sum(cT)
                cell_dict[cell]['c'] = np.sum(g)
                cell_dict[cell]['g'] = np.sum(c)
    ragged_h5_local.close()
    q.put(cell_dict)
    return None
def make_write_function(h5file):
    def write_hdf5(q):
        full_cell_dict = {}
        while True:
            cell_d = q.get()
            if cell_d is None: break
            for cell,d in cell_d.items():
                if cell in full_cell_dict:
                    for k,v in d.items():
                        if k in ['gene','sc','tc']:
                            full_cell_dict[cell][k].extend(v)
                        else:
                            full_cell_dict[cell][k] += v
                else:
                    full_cell_dict[cell] = {}
                    for k,v in d.items():
                        full_cell_dict[cell][k] = v

            del cell_d
            q.task_done()
        hdf5_file = h5py.File(h5file, 'a')
        grp = hdf5_file.require_group('cells')
        for cell, d in full_cell_dict.items():
            cell_grp = grp.require_group(cell)
            for tag, value in d.items():
                if tag == 'gene':
                    cell_grp.create_dataset(tag, data=np.array(value, dtype='S'))
                elif tag in ['sc','tc']:
                    dset = cell_grp.create_dataset(tag, (len(value),),dtype=h5py.vlen_dtype(np.dtype('int32')))
                    for i,arr in enumerate(value):
                        dset[i] = arr
                else:
                    cell_grp.create_dataset(tag,data=value)
        hdf5_file.close()
        q.task_done()
    return write_hdf5
def estim_Pc(cell_id, h5file):
    fd = h5py.File(h5file, 'r')
    grp = fd['cells/{}'.format(cell_id)]
    A = sparse.dok_matrix((500,500), dtype=np.int32)
    for k,n in zip(np.concatenate(grp['sc'][:]),np.concatenate(grp['tc'][:])):
        if k < 500 and n < 500:
            A[k,n] += 1
    A = A.toarray()
    p_e = (grp['gA'][()]/grp['g'][()] + grp['cT'][()]/grp['c'][()])/2
    fd.close()
    return estimateP_c(A,p_e,cell_id)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract molecule information and store in an hdf5 file')
    parser.add_argument('-h5','--hdf5',metavar='input', type=str, help='.h5 file to process')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--skip_copy',action='store_true')
    args = parser.parse_args()
    h5file = args.hdf5
    threads = args.threads
    skip_copy = args.skip_copy
    if not skip_copy:
        print('Copying over gene information to cell information')
        fd = h5py.File(h5file,'r')
        gene_list = list(fd['genes'].keys())
        fd.close()
        m = Manager()
        q = m.JoinableQueue()
        p = Process(target=make_write_function(h5file), args=(q,))
        p.start()
        res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(get_cell_dicts)(genes,h5file,q) for genes in np.array_split(gene_list,threads))
        q.put(None)
        print('Done, writing results')
        p.join()
    print('Estimating p_c and p_e')
    fd = h5py.File(h5file,'r')
    cell_list = list(fd['cells'].keys())
    fd.close()
    params = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(estim_Pc)(cell_id, h5file) for cell_id in cell_list)
    res_dict = {cell_id: (p_c,p_e) for p_c,p_e,cell_id in params}
    print('Done, writing results')
    fd = h5py.File(h5file,'a')
    for cell_id, cell_grp in fd['cells'].items():
        cell_grp.attrs['p_c'] = res_dict[cell_id][0]
        cell_grp.attrs['p_e'] = res_dict[cell_id][1]
    fd.close()

