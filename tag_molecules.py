import pysam
import numpy as np
import pandas as pd
import portion as P
import itertools
import pygtrie
import vcf
import gc
import argparse
import seaborn as sns
from collections import Counter
from pyfaidx import Fasta
from scipy.stats import binom
import os
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
import h5py
import time

# taken from https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]
        
def compare(read_seq, ref_seq, is_reverse, start):
    specific_conversions = {}
    total_content = {'a' : ref_seq.count('A'), 'c' : ref_seq.count('C'), 'g' : ref_seq.count('G'), 't' : ref_seq.count('T')}
    specific_conversions[('c', 'A')] = 0
    specific_conversions[('g', 'A')] = 0
    specific_conversions[('t', 'A')] = 0
    specific_conversions[('a', 'C')] = 0
    specific_conversions[('g', 'C')] = 0
    specific_conversions[('t', 'C')] = 0
    specific_conversions[('a', 'G')] = 0
    specific_conversions[('c', 'G')] = 0
    specific_conversions[('t', 'G')] = 0
    specific_conversions[('a', 'T')] = 0
    specific_conversions[('c', 'T')] = 0
    specific_conversions[('g', 'T')] = 0
    specific_conversions[('a', 'N')] = 0
    specific_conversions[('c', 'N')] = 0
    specific_conversions[('g', 'N')] = 0
    specific_conversions[('t', 'N')] = 0
    specific_conversions[('n', 'A')] = 0
    specific_conversions[('n', 'C')] = 0
    specific_conversions[('n', 'T')] = 0
    specific_conversions[('n', 'G')] = 0
    conv_loc = []
    mut_loc = {}
    for i, (read_base, ref_base) in enumerate(zip(read_seq, ref_seq)):
        if read_base == ref_base:
            continue
        else:
            k = (ref_base.lower(),read_base)
            specific_conversions[k] += 1
            if is_reverse:
                if ''.join(k) == 'aG':
                    conv_loc.append(start + i)
            else:
                if ''.join(k) == 'tC':
                    conv_loc.append(start + i)
            if ref_base.lower() == 'n':
                mut_loc[start+i] = read_base
    return Counter(total_content), Counter(specific_conversions), conv_loc, mut_loc
def get_allele(mol, mut_trie, comparison_dict, strand):
    ref_count = 0
    alt_count = 0
    for p in mut_trie[mol].keys():
            try:def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
                comparison_dict[p]
            except KeyError:
                continue
            if mut_trie[mol][p] == comparison_dict[p][0]:
                ref_count += 1
            elif mut_trie[mol][p] == comparison_dict[p][1]:
                alt_count += 1
            else:
                if strand == '+' and mut_trie[mol][p] == 'C':
                    ref_count += int(comparison_dict[p][0] == 'T')
                    alt_count += int(comparison_dict[p][1] == 'T')
                elif mut_trie[mol][p] == 'G':
                    ref_count += int(comparison_dict[p][0] == 'A')
                    alt_count += int(comparison_dict[p][1] == 'A')
    return mol, ref_count, alt_count
def get_alleles(mols,mut_trie,comparison_dict,strand):
    for mol in mols:
        mol, ref_count, alt_count = get_allele(mol, mut_trie, comparison_dict, strand)
        if ref_count == 0 and alt_count == 0:
            yield np.nan
        elif ref_count > 0 and alt_count == 0:
            yield True
        elif alt_count > 0 and ref_count == 0:
            yield False
        else:
            yield np.nan

def compare_to_allele(mutation_allele, mol,mut_trie,comparison_dict,strand):
    mol, ref_count, alt_count = get_allele(mol,mut_trie,comparison_dict,strand)
    if ref_count == 0 and alt_count == 0:
        return mol, True
    elif ref_count > 0 and alt_count == 0:
        is_ref = True
    elif alt_count > 0 and ref_count == 0:
        is_ref = False
    else:
        is_ref = ref_count > alt_count
    return mol,is_ref == mutation_allele

def get_coverage(p, intervals_list):
    return np.sum(np.array([p in interval for interval in intervals_list]))

def get_num_cells(loc):
    return len(set([k.partition(':')[0] for k in loc.index[loc]]))

def createTag(d):
    return ''.join([''.join(key) + str(d[key]) + ';' for key in d.keys()])[:-1]

def get_mols_w_mut(loc,ref,mut_trie,comparison_dict,strand, locs_df):
    mols = []
    for r in locs_df[locs_df[loc]].apply(lambda row: compare_to_allele(ref,row.name,mut_trie,comparison_dict,strand), axis=1):
        if r[1]:
            mols.append(r[0])
    return mols

def get_tries(bam,fasta_file, g_dict):
    fasta_ref = Fasta(fasta_file)
    read_trie = pygtrie.StringTrie(separator=':')
    conv_trie = pygtrie.StringTrie(separator=':')
    mut_trie = pygtrie.StringTrie(separator=':')
    t_trie = pygtrie.StringTrie(separator=':')
    s_trie = pygtrie.StringTrie(separator=':')
    strand = pygtrie.StringTrie(separator=':')
    for read in bam.fetch(g_dict['seqid'],g_dict['start'],g_dict['end']):
        if read.get_tag('XT') != g_dict['gene_id']:
            continue

        read_trie[read.query_name] = read
        offset = 0
        tc_list = []
        sc_list = []
        cl_list = []
        mut_list = []
        intervals = intervals_extract(read.get_reference_positions())
        for inter in intervals:
            ref_seq = fasta_ref[read.reference_name][inter[0]:inter[1]].seq
            read_seq = read.query_alignment_sequence[offset:offset + (inter[1]-inter[0])]
            offset += inter[1]-inter[0] + 1
            tc, sc, cl, ml = compare(read_seq, ref_seq, read.is_reverse, inter[0])
            tc_list.append(tc)
            sc_list.append(sc)
            cl_list.append(cl)
            mut_list.append(ml)
        strand[read.query_name] = read.is_reverse
        t = Counter()
        [t.update(i) for i in tc_list]
        s = Counter()
        [s.update(i) for i in sc_list]
        full_cl_list = [item for sublist in cl_list for item in sublist]
        full_mut_list = {item[0]:item[1] for sublist in mut_list for item in sublist.items()}
        t_trie[read.query_name] = t
        s_trie[read.query_name] = s
        conv_trie[read.query_name] = (full_cl_list, [i for i in intervals_extract(read.get_reference_positions())])
        mut_trie[read.query_name] = full_mut_list
    fasta_ref.close()
    return read_trie,mut_trie,conv_trie,t_trie, s_trie,strand

def find_mutations(conv_trie,mut_trie,g_dict, vcf_reader):
    locs_list = []
    intervals_list = []
    for cell,convs in sorted(conv_trie.items()):
        if len(convs[0]) > 0:
            s = pd.Series(np.repeat(True,len(convs[0])), index=convs[0])
            s.name = cell
            locs_list.append(s)
        intervals_list.append(P.from_data([[True,p[0],p[1],True] for p in convs[1]]))

    locs_df = pd.DataFrame(locs_list).fillna(False)
    pos_coverage = locs_df.columns.map(lambda x: get_coverage(int(x), intervals_list))
    cov_series = pd.Series(pos_coverage.values,index=locs_df.columns)

    stats_df = pd.DataFrame([locs_df.sum(),cov_series,(locs_df.sum()/cov_series),locs_df.apply(get_num_cells)]).T
    stats_df.columns = ['conversions','coverage','fraction','num_cells']

    if stats_df.empty:
        return False, None, None

    comparison_dict = {}
    if vcf_file is not None:
        vcf_reader = vcf.Reader(filename=vcf_file)
        for i in vcf_reader.fetch(g_dict['seqid'],g_dict['start'],g_dict['end']):
            comparison_dict[i.POS-1] = (i.REF, i.ALT[0])
        del vcf_reader
    p_median = stats_df['fraction'].median()
    stats_df['pval'] = stats_df.apply(lambda row: binom.sf(row['conversions'],row['coverage'], p_median), axis=1)
    if vcf_file is not None:
        frac_ref = stats_df.apply(lambda row: np.array(list(get_alleles(locs_df.index[locs_df[row.name]].values,mut_trie,comparison_dict,g_dict['strand']))), axis=1).apply(lambda x: np.mean(x[~np.isnan(x)]) if np.mean(np.isnan(x)) != 1 else np.nan)
        frac_ref.name = 'frac_ref'
        stats_df = stats_df.join(frac_ref)
        sig_df = stats_df[(stats_df['num_cells']>4).mul((stats_df['frac_ref'] - 0.5).abs() > 0.3).mul(stats_df['pval'] < 0.05/stats_df.shape[0])].sort_values('pval')
        sig_df['ref'] = sig_df['frac_ref'] > 0.5
    sig_df = stats_df[(stats_df['num_cells']>4).mul(stats_df['pval'] < 0.05/stats_df.shape[0])].sort_values('pval')
    sig_df['ref'] = True
    sig_df['pos'] = sig_df.index
    sig_df['chrom'] = g_dict['seqid']
    sig_df['gene'] = g_dict['gene_id']
    mols_series = sig_df.apply(lambda row: get_mols_w_mut(row.name,row['ref'],mut_trie,comparison_dict,g_dict['strand'], locs_df), axis=1)
    return True, sig_df, mols_series

def count_conversions(bfile,fasta_file,vcf_file,g_dict,q):
    bam = pysam.AlignmentFile(bfile)
    read_trie,mut_trie,conv_trie,t_trie, s_trie,strand_trie = get_tries(bam,fasta_file,g_dict)
    good, sig_df, mols_series = find_mutations(conv_trie,mut_trie,g_dict, vcf_file)
    if not good:
        return g_dict['gene_id']
    for pos, mols in mols_series.items():
        for mol in mols:
            if strand_trie[mol]:
                key_tuple = ('a', 'G')
            else:
                key_tuple = ('t', 'C')
            conv_trie[mol] = (list(set(conv_trie[mol][0]) - set([pos])), conv_trie[mol][1])
            s_trie[mol][key_tuple] -= 1
            t_trie[mol][key_tuple[0]] -= 1
            t_trie[mol][key_tuple[1].lower()] += 1
    content_list = []
    append_content = content_list.append
    for name,mol in read_trie.iteritems():
        mol.set_tag('SC',createTag(s_trie[mol.query_name]),'Z')
        mol.set_tag('TC',createTag(t_trie[mol.query_name]),'Z')
        if len(conv_trie[mol.query_name][0]) > 0:
            mol.set_tag('CL',conv_trie[mol.query_name][0])
        content_dict = {'mol_string': mol.to_string()}
        append_content(content_dict)
    # chunking to avoid trying to send too much data
    if len(content_list) > 10000:
        for c_list in chunks(content_list, 10000):
            q.put((True,c_list))
    else:
        q.put((True,content_list))
    bam.close() 
    if sig_df.empty:
        return g_dict['gene_id']
    cell_series = mols_series.apply(lambda row: list(set([m.partition(':')[0] for m in row])))
    cell_series.name = 'cells'
    sig_dict = sig_df.join(cell_series).T.to_dict()
    q.put((False, sig_dict))
    return g_dict['gene_id']

def create_h5_function(infile, outfile, h5outfile, version):
    bamfile = pysam.AlignmentFile(infile, 'rb')
    h = bamfile.header
    def write_new_reads_and_pc_data(q):
        new_bamfile = pysam.AlignmentFile(outfile, mode='wb',template=bamfile)
        f = h5py.File(h5outfile, 'a', libver='latest')
        while True:
            is_mol,content_list = q.get()
            if is_mol == None: break
            if is_mol:
                for content_dict in content_list:
                    new_bamfile.write(pysam.AlignedRead.fromstring(content_dict['mol_string'],h))
            else:
                 for pos,d in content_list.items():
                    g = f.require_group('mutations')
                    gg = g.require_group('{chrom}/{gene}/{pos}/'.format(chrom=d['chrom'],gene=d['gene'], pos=d['pos']))
                    dset = gg.create_dataset('conversions', data=np.int_(d['conversions'])) 
                    dset = gg.create_dataset('coverage', data=np.int_(d['coverage']))
                    dset = gg.create_dataset('fraction', data=d['fraction'])
                    dset = gg.create_dataset('pval', data=d['pval'])
                    dset = gg.create_dataset('ref', data=d['ref'])
                    dset = gg.create_dataset('cells', data=np.array(d['cells'],dtype='S'))
                    dset = gg.create_dataset('gene', data=d['gene'])
            q.task_done()
        f.close()
        q.task_done()
        return None
    return write_new_reads_and_pc_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tag molecules for NASC-seq2')
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output',type=str, help='Output .bam file')
    parser.add_argument('-g','--gtf',type=str, help='gtf file with gene information')
    parser.add_argument('-f','--fasta',type=str, help='Fasta file for your reference genome')
    parser.add_argument('-v','--vcf',type=str,default=None, help='vcf file with genotype information')
    parser.add_argument('-mut','--mutations',type=str, help='Output hdf5 file with mutation information')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--contig', default=None, metavar='contig', type=str, help='Restrict stitching to contig')
    args = parser.parse_args()
#    vcf_file = '/mnt/storage1/home/antonl/meta/vcf/CAST.SNPs.validated.vcf.gz'
    vcf_file = args.vcf
#    gtffile = '/mnt/davidson/hendgert/recovered/resources/gtf/mouse/Mus_musculus.GRCm38.91.chr.clean.gtf'
    gtffile = args.gtf
#    bfile = '/mnt/storage1/home/antonl/projects/NASC-seq/testing/mESC_NASCseq_EXP-20-CB7751.filtered.Aligned.GeneTagged.UBcorrected.sorted.stitched.fixed.sorted.bam'
    bfile = args.input
#    fasta_file = '/mnt/davidson/hendgert/recovered/resources/genomes/Mouse_CAST_Nmasked/mm10_Nmasked.fa'
    fasta_file = args.fasta
#    outfile = 'mESC_NASCseq_EXP-20-CB7751.sam'
    outfile = args.output
#    h5outfile = 'mESC_NASCseq_EXP-20-CB7751.h5'
    h5outfile = args.mutations
   
#    contig = None
    contig = args.contig
    threads = int(args.threads)
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
    gene_df = pd.DataFrame(gene_list)
    gene_df.index = gene_df['gene_id']
    del gene_list
    m = Manager()
    q = m.JoinableQueue()
    p = Process(target=create_h5_function(infile=bfile, outfile=outfile, h5outfile=h5outfile, version='0.1'), args=(q,))
    p.start()

    print('Counting conversions in molecules for {}'.format(bfile))

    start = time.time()
    params = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(count_conversions)(bfile,fasta_file,vcf_file, gene,q) for g,gene in gene_df.iterrows())
    q.put((None,None))
    p.join()
    end = time.time()
    print(end-start)
    print('done')
