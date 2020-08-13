import pysam
import argparse
import h5py
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
import time
import pygtrie

def parseTCTag(read):
    tag = read.get_tag('TC')
    splittag = tag.split(';')
    total_content = {}
    for c in splittag:
        total_content[c[0]] = np.int_(c[1:])
    return total_content

def parseSCTag(read):
    tag = read.get_tag('SC')
    splittag = tag.split(';')
    specific_conversions = {}
    for c in splittag:
        specific_conversions[(c[0], c[1])] = np.int_(c[2:])
    return specific_conversions

def get_tags(filename_bam, g_dict, q):
    content_dict = {}
    bam = pysam.AlignmentFile(filename_bam, 'r')
    for mol in bam.fetch(g_dict['seqid'],g_dict['start'],g_dict['end']):
        cell =  mol.get_tag('BC')
        gene = mol.get_tag('XT')
        umi = mol.get_tag('UB')
        if g_dict['gene_id'] != mol.get_tag('XT'):
            continue
        content = {'umi': umi,
                'sc': {''.join(k):v for k,v in parseSCTag(mol).items()},
                'tc': {''.join(k):v for k,v in parseTCTag(mol).items()},
                'nr': mol.get_tag('NR'),
                'ir': mol.get_tag('IR'),
                'er': mol.get_tag('ER'),
                'start': mol.reference_start,
                'end': mol.reference_end,
                'rl': mol.query_length
                }
        if mol.has_tag('CL'):
            content['cl'] = mol.get_tag('CL')
            content['cl_name'] = '{}/{}'.format(cell,umi)
        if cell in content_dict:
            content_dict[cell].append(content)
        else:
            content_dict[cell] = [content]
    bam.close()
    if len(content_dict) == 0:
        return g_dict['gene_id']
    nicer_content_dict = {'cell': [], 'umi': [], 'sc': {'cA': [],
   'gA': [],
   'tA': [],
   'aC': [],
   'gC': [],
   'tC': [],
   'aG': [],
   'cG': [],
   'tG': [],
   'aT': [],
   'cT': [],
   'gT': [],
   'aN': [],
   'cN': [],
   'gN': [],
   'tN': [],
   'nA': [],
   'nC': [],
   'nT': [],
   'nG': []},
                      'tc': {'a':[], 'c': [], 'g':[], 't': []},
                     'nr': [], 'ir':[], 'er':[], 'start': [],
                     'end': [], 'rl':[],'cl':[], 'cl_name':[]}
    for cell,cell_dicts in content_dict.items():
        nicer_content_dict['cell'].append(cell)
        temp_dict = {}
        for d in cell_dicts:

            for tag, v in d.items():
                if tag not in temp_dict:
                    if tag in ['sc', 'tc']:
                        temp_dict[tag] = {}
                        for k,v2 in v.items():
                            temp_dict[tag][k] = [v2]
                    else:
                        temp_dict[tag] = [v]
                else:
                    if tag in ['sc', 'tc']:
                        for k,v2 in v.items():
                            temp_dict[tag][k].append(v2)
                    else:
                        temp_dict[tag].append(v)

        for k,v in temp_dict.items():
            if k in ['sc', 'tc']:
                for k2, v2 in v.items():
                    nicer_content_dict[k][k2].append(v2)
            elif k in ['cl','cl_name']:
                nicer_content_dict[k].extend(v)
            else:
                if k == 'umi':
                    nicer_content_dict[k].append(v)
                else:
                    nicer_content_dict[k].append(v)
    q.put((g_dict, nicer_content_dict))
    return g_dict['gene_id']

def create_h5_function(h5outfile):
    def write_new_reads_and_pc_data(q):
        f = h5py.File(h5outfile, 'a', libver='latest')
        grp = f.create_group('genes')
        while True:
            gene_dict,content_dict = q.get()
            if content_dict == None: break
            gene_grp = grp.create_group(gene_dict['gene_id'])
            gene_grp.attrs['start'] = gene_dict['start']
            gene_grp.attrs['end'] = gene_dict['end']
            gene_grp.attrs['strand'] = gene_dict['strand']
            gene_grp.attrs['chrom'] = gene_dict['seqid']
            for tag, value in content_dict.items():
                if tag == 'umi':
                    dt = h5py.special_dtype(vlen=str)
                else:
                    dt = h5py.vlen_dtype(np.dtype('int32'))
                if tag in ['sc', 'tc']:
                    tag_grp = gene_grp.create_group(tag)
                    for tag2, v2 in value.items():
                        dataset = tag_grp.create_dataset(tag2,(len(v2),),dtype=dt)
                        for i,arr in enumerate(v2):
                            dataset[i] = arr
                elif tag in ['cell','cl_name']:
                    dataset = gene_grp.create_dataset(tag, data=np.array(value, dtype='S'))
                elif tag == 'umi':
                    dataset = gene_grp.create_dataset(tag, data=np.array(value, dtype=object), dtype=dt)
                else:
                    dataset = gene_grp.create_dataset(tag, (len(value),), dtype=dt)
                    for i, arr in enumerate(value):
                        dataset[i] = arr
            q.task_done()
        f.close()
        q.task_done()
        return None
    return write_new_reads_and_pc_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract molecule information and store in an hdf5 file')
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .h5 file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--contig', default=None, metavar='contig', type=str, help='Restrict extraction to contig')
    args = parser.parse_args()
    #filename_bam = 'mESC_NASCseq_EXP-20-CB7751.sorted.bam'
    filename_bam = args.input
    #filename_h5 = 'mESC_NASCseq_EXP-20-CB7751.sorted.h5'
    filename_h5 = args.output
    #gtffile = '/mnt/davidson/hendgert/recovered/resources/gtf/mouse/Mus_musculus.GRCm38.91.chr.clean.gtf'
    gtffile = args.gtf
    #threads = 300
    threads = args.threads
    #contig = None
    contig = args.contig
    gene_list = []
    with open(gtffile, 'r') as f:
        for line in f:
            l = line.split('\t')
            if len(l) < 8:
                continue
            if l[2] == 'gene':
                if contig is not None:
                    if l[0] == contig:
                        gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                    else:
                        continue
                else:
                    gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
    gene_dict = {g['gene_id']: g for g in gene_list}
    print('Saving tags for {}'.format(filename_bam))

    start = time.time()
    m = Manager()
    q = m.JoinableQueue()
    p = Process(target=create_h5_function(h5outfile=filename_h5), args=(q,))
    p.start()
    params = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(get_tags)(filename_bam, g_dict,q) for g,g_dict in gene_dict.items())
    q.put((None,None))
    p.join()
    end = time.time()
    print(end-start)
    print('done') 
