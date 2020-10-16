import pysam
import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
from joblib import Parallel, delayed
from time import time


def add_umi_tags(infile, outfile, threads):
    bam_in = pysam.AlignmentFile(infile, 'rb', threads=int(threads / 2))
    bam_out = pysam.AlignmentFile(outfile, 'wb', template=bam_in,threads=int(threads / 2))
    condition_gene_dict = {}
    for read in bam_in.fetch():
        # Do we need a try statement here? What about phi-X reads?
        condition = read.get_tag('BC')
        try:
            gene = read.get_tag('GE')
        except KeyError:
            continue
        try:
            condition_gene_dict[condition][gene] += 1
        except KeyError:
            try:
                condition_gene_dict[condition].update({gene:1})
            except KeyError:
                condition_gene_dict.update({condition:{gene:1}})
        umi = condition_gene_dict[condition][gene]
        read.set_tag('UB', value = str(umi), value_type = 'Z')
        bam_out.write(read)

if __name__ == "__main__":
    tic = time()
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-p', '--n_jobs', type = int)
    o = parser.parse_args()

    infile = o.infile
    outfile = o.outfile
    n_jobs = o.n_jobs

    add_umi_tags(infile, outfile, n_jobs)

    toc = time()

    print(toc-tic)
