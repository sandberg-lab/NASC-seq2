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
    f = h5py.File(h5outfile + gene + '.h5', 'a')
    print(gene, len(content_list))
    for content_dict in content_list:
        name = content_dict['name']
        if 'cells/{name}'.format(name=name) in f:
            continue
        f.create_dataset('cells/{name}/sc/cA'.format(name=name), data=content_dict['sc'][('c','A')])
        f.create_dataset('cells/{name}/sc/gA'.format(name=name), data=content_dict['sc'][('g','A')])
        f.create_dataset('cells/{name}/sc/tA'.format(name=name), data=content_dict['sc'][('t','A')])
        f.create_dataset('cells/{name}/sc/gC'.format(name=name), data=content_dict['sc'][('g','C')])
        f.create_dataset('cells/{name}/sc/tC'.format(name=name), data=content_dict['sc'][('t','C')])
        f.create_dataset('cells/{name}/sc/aC'.format(name=name), data=content_dict['sc'][('a','C')])
        f.create_dataset('cells/{name}/sc/aG'.format(name=name), data=content_dict['sc'][('a','G')])
        f.create_dataset('cells/{name}/sc/cG'.format(name=name), data=content_dict['sc'][('c','G')])
        f.create_dataset('cells/{name}/sc/tG'.format(name=name), data=content_dict['sc'][('t','G')])
        f.create_dataset('cells/{name}/sc/aT'.format(name=name), data=content_dict['sc'][('a','T')])
        f.create_dataset('cells/{name}/sc/cT'.format(name=name), data=content_dict['sc'][('c','T')])
        f.create_dataset('cells/{name}/sc/gT'.format(name=name), data=content_dict['sc'][('g','T')])
        f.create_dataset('cells/{name}/sc/aN'.format(name=name), data=content_dict['sc'][('a','N')])
        f.create_dataset('cells/{name}/sc/cN'.format(name=name), data=content_dict['sc'][('c','N')])
        f.create_dataset('cells/{name}/sc/gN'.format(name=name), data=content_dict['sc'][('g','N')])
        f.create_dataset('cells/{name}/sc/tN'.format(name=name), data=content_dict['sc'][('t','N')])
        f.create_dataset('cells/{name}/sc/nA'.format(name=name), data=content_dict['sc'][('n','A')])
        f.create_dataset('cells/{name}/sc/nC'.format(name=name), data=content_dict['sc'][('n','C')])
        f.create_dataset('cells/{name}/sc/nT'.format(name=name), data=content_dict['sc'][('n','T')])
        f.create_dataset('cells/{name}/sc/nG'.format(name=name), data=content_dict['sc'][('n','G')])
        for b in bases:
            f.create_dataset('cells/{name}/tc/{base}'.format(name=name,base=b), data=content_dict['tc'][b])
        if 'cl' in content_dict:
            f.create_dataset('cells/{name}/cl'.format(name=name), data=content_dict['cl'])
        f.create_dataset('cells/{name}/nr'.format(name=name), data=content_dict['nr'])
        f.create_dataset('cells/{name}/ir'.format(name=name), data=content_dict['ir'])
        f.create_dataset('cells/{name}/er'.format(name=name), data=content_dict['er'])
        f.create_dataset('cells/{name}/start'.format(name=name), data=content_dict['start'])
        f.create_dataset('cells/{name}/end'.format(name=name), data=content_dict['end'])
        f.create_dataset('cells/{name}/rl'.format(name=name), data=content_dict['rl'])
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
        if g_dict['gene_id'] != mol.get_tag('XT'):
            continue
        name = mol.query_name
        name = name.replace(':','/')
        content_dict = {'name': name,
                'sc': parseSCTag(mol),
                'tc': parseTCTag(mol),
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
        write_tags(content_list, g_dict['gene_id'])
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
