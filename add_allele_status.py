import pygtrie
import argparse
import tensorflow as tf
import pysam
import h5py
import vcf
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
from multiprocessing import Process, Manager
from joblib import delayed,Parallel

def get_allele(hits, comparison_dict, strand, new_mol):
    ref_count = 0
    alt_count = 0
    for p in hits.keys():
            try:
                comparison_dict[p]
            except KeyError:
                continue
            if hits[p] == comparison_dict[p][0]:
                ref_count += 1
            elif hits[p] == comparison_dict[p][1]:
                alt_count += 1
            else:
                if strand == '+' and hits[p] == 'C' and new_mol:
                    ref_count += int(comparison_dict[p][0] == 'T')
                    alt_count += int(comparison_dict[p][1] == 'T')
                elif strand =='-' and hits[p] == 'G' and new_mol:
                    ref_count += int(comparison_dict[p][0] == 'A')
                    alt_count += int(comparison_dict[p][1] == 'A')
    return ref_count, alt_count

def get_alleles_per_gene_no_mutation_file(gene, h5_filename, bamfile_stitched, vcf_file, q):
    h5_file = h5py.File(h5_filename,'r')
    stitcher_bam = pysam.AlignmentFile(bamfile_stitched, 'r')

    new_trie = pygtrie.StringTrie()
    example_gene_h5 = h5_file['genes/{}'.format(gene)]
    s = example_gene_h5.attrs['strand']
    cell_list = [cell.decode('utf-8') for cell in example_gene_h5['cell'][:]]
    try:
        cell_dict = {cell: eval(umis) for cell,umis in zip(cell_list,example_gene_h5['umi'][:])}
    except:
        cell_dict = {cell: umis for cell,umis in zip(cell_list,example_gene_h5['umi'][:])}
    htest_dict = {cell: htests for cell, htests in zip(cell_list,example_gene_h5['htest'][:])}
    for cell, l in cell_dict.items():
        htests = htest_dict[cell]
        for umi,new in zip(l,htests):
            new_trie['{}/{}'.format(cell,umi)] = new

    vcf_reader = vcf.Reader(filename=vcf_file)

    comparison_dict = {}
    for i in vcf_reader.fetch(example_gene_h5.attrs['chrom'], example_gene_h5.attrs['start'],example_gene_h5.attrs['end']):
            comparison_dict[i.POS-1] = (i.REF, i.ALT[0])

    allele_position_set = set(comparison_dict.keys())
    
    allele_dict = {}
    for mol in stitcher_bam.fetch(example_gene_h5.attrs['chrom'], example_gene_h5.attrs['start'],example_gene_h5.attrs['end']):
        cell = mol.get_tag('BC')
        umi = mol.get_tag('UB')
        if gene != mol.get_tag('XT'):
                continue
        reference_positions = mol.get_reference_positions()
        read_sequence = mol.query_alignment_sequence
        hits = {x:read_sequence[i]  for i,x in enumerate(reference_positions) if x in allele_position_set}

        if len(hits) > 0:
            try:
                a = get_allele(hits,comparison_dict,s,new_trie['{}/{}'.format(cell,umi)])
                count = a
            except KeyError:
                count = (0, 0)
        else:
            count = (0,0)
        if count[0] > 0 and count[1] == 0:
            res = 1.0
        elif count[1] > 0 and count[0] == 0:
            res = 0.0
        else:
            res = np.nan
        if cell in allele_dict:
            allele_dict[cell].append(res)
        else:
            allele_dict[cell] = [res]
    
    new_ragged = tf.cast(tf.ragged.constant(example_gene_h5['htest'][:]), dtype=float)
    allele_list = [allele_dict[cell] for cell in cell_list]
    allele_ragged = tf.ragged.constant(allele_list)
    ref_new = allele_ragged*new_ragged
    alt_new = (1-allele_ragged)*new_ragged
    
    res_dict = {'allele_call': allele_ragged, 'ref_new': ref_new, 'alt_new': alt_new}
    f_w = h5py.File('{}_tmp.h5'.format(gene), 'w', libver='latest')
    gene_grp = f_w.create_group(gene)
    dt = h5py.vlen_dtype(np.dtype('int32'))
    for tag, value in res_dict.items():
        dset = gene_grp.create_dataset(tag,  (value.shape[0],), dtype=dt)
        for i, arr in enumerate(value):
            dset[i] = arr
    f_w.close()
    h5_file.close()
    stitcher_bam.close()
    q.put(gene)
    return gene
def create_h5_function(h5outfile):
    def write_h5_file(q):
        f = h5py.File(h5outfile, 'a', libver='latest')
        grp = f.create_group('genes')
        while True:
            gene_name = q.get()
            if gene_name == None: break
            f_gene = h5py.File('{}_tmp.h5'.format(gene_name), 'r', libver='latest')
            f_gene.copy(gene_name, grp)
            f_gene.close()
            os.system("rm {}_tmp.h5".format(gene_name))
            q.task_done()
        f.close()
        q.task_done()
        return None
    return write_h5_file

def main():
    parser = argparse.ArgumentParser(description='Add allele status per molecule')
    parser.add_argument('-h5_in','--hdf5_in',metavar='input', type=str, help='.h5 file to process')
    parser.add_argument('-b', '--bamfile', metavar='bam', type=str, default=1, help='Stitched bam file')
    parser.add_argument('-v', '--vcf_file', metavar='vcf', type=str, default=1, help='vcf file')
    parser.add_argument('-h5_out','--hdf5_out',metavar='output', type=str, help='.h5 file to process')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')


    args = parser.parse_args()
    
    h5_filename = args.hdf5_in
    bamfile_stitched = args.bamfile
    vcf_file = args.vcf_file
    h5_out = args.hdf5_out
    threads = args.threads

    h5_file = h5py.File(h5_filename, 'r')
    gene_list = list(h5_file['genes'].keys())
    h5_file.close()

    m = Manager()
    q = m.JoinableQueue()
    p = Process(target=create_h5_function(h5outfile=h5_out), args=(q,))
    print('Adding allele status to each molecule, {} threads'.format(threads))
    p.start()
    params = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(get_alleles_per_gene_no_mutation_file)(gene, h5_filename, bamfile_stitched, vcf_file, q) for gene in gene_list)
    q.put(None)
    p.join()
    print('Done.')

if __name__ == '__main__':
    main()
