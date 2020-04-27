import pysam
import h5py
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
import time
import pygtrie
def write_tags(content_list, gene):
    bases = ['a','t','c','g']
    f = h5py.File(h5outfile + gene + '_new.h5', 'a')
    print(gene, len(content_list))
    for content_dict in content_list:
        name = content_dict['name']
        if '/{name}'.format(name=name) in f:
            continue
        grp = f.create_group('/{name}'.format(name=name))
        grp.create_dataset('sc'.format(name=name), data=content_dict['sc'])
        grp.create_dataset('tc'.format(name=name), data=content_dict['tc'])
        if 'cl' in content_dict:
            grp.create_dataset('cl'.format(name=name), data=content_dict['cl'])
        
        grp.attrs['nl'] = content_dict['nr']
        grp.attrs['ir'] = content_dict['ir']
        grp.attrs['er'] = content_dict['er']
        grp.attrs['start'] = content_dict['start']
        grp.attrs['end'] = content_dict['end']
        grp.attrs['rl'] = content_dict['rl']
    f.close()
    return None

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


def get_tags(filename_bam, g_dict):
    content_list = []
    content_append = content_list.append
    bam = pysam.AlignmentFile(filename_bam, 'r')
    for mol in bam.fetch(g_dict['seqid'],g_dict['start'],g_dict['end']):
        cell =  mol.get_tag('BC')
        gene = mol.get_tag('XT')
        umi = mol.get_tag('UB')
        if g_dict['gene_id'] != mol.get_tag('XT'):
            continue
        name = '{}/{}'.format(cell,umi)
        content_dict = {'name': name,
                'sc': list(parseSCTag(mol).values()),
                'tc': list(parseTCTag(mol).values()),
                'nr': mol.get_tag('NR'),
                'ir': mol.get_tag('IR'),
                'er': mol.get_tag('ER'),
                'start': mol.reference_start,
                'end': mol.reference_end,
                'rl': mol.query_length
                }
        if mol.has_tag('CL'):
            content_dict['cl'] = mol.get_tag('CL')
        content_append(content_dict)
    if len(content_list) > 0:
        write_tags_new(content_list, g_dict['gene_id'])
    bam.close()
    return g_dict['gene_id']

if __name__ == '__main__':
    filename_bam = 'mESC_NASCseq_EXP-20-CB7751.sorted.bam'
    filename_h5 = 'mESC_NASCseq_EXP-20-CB7751.sorted.h5'
    gtffile = '/mnt/davidson/hendgert/recovered/resources/gtf/mouse/Mus_musculus.GRCm38.91.chr.clean.gtf'
    h5outfile = 'h5_files/mESC_NASCseq_EXP-20-CB7751.'
    threads = 300
    contig = None
    gene_list = []
    with open(gtffile, 'r') as f:
        for line in f:
            l = line.split('\t')
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
    params = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(get_tags)(filename_bam, gene) for g,gene in gene_dict.items())
    end = time.time()
    print(end-start)
    print('done') 
