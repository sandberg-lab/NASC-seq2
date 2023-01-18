import matplotlib, argparse
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed

from scipy.stats import spearmanr

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

parser = argparse.ArgumentParser()
parser.add_argument('--load', action='store_true')
o = parser.parse_args()

def sns_styleset():
    sns.set_context('paper')
    sns.set_style('ticks')
    matplotlib.rcParams['axes.linewidth']    = .75
    matplotlib.rcParams['xtick.major.width'] = .75
    matplotlib.rcParams['ytick.major.width'] = .75
    matplotlib.rcParams['xtick.major.size'] = 3
    matplotlib.rcParams['ytick.major.size'] = 3
    matplotlib.rcParams['xtick.minor.size'] = 2
    matplotlib.rcParams['ytick.minor.size'] = 2
    matplotlib.rcParams['font.size']       = 7
    matplotlib.rcParams['axes.titlesize']  = 7
    matplotlib.rcParams['axes.labelsize']  = 7
    matplotlib.rcParams['legend.fontsize'] = 7
    matplotlib.rcParams['xtick.labelsize'] = 7
    matplotlib.rcParams['ytick.labelsize'] = 7

sns_styleset()


only_x = False
same_strand = False
discretize = False
merge_alleles = False
cellt = 500
discretize_t = 3

def compute_correlation(x,y):
    if discretize:
        x = [val > discretize_t and 1 or 0 for val in x]
        y = [val > discretize_t and 1 or 0 for val in y]
        res = spearmanr(x,y)[0]
    else:
        res = spearmanr(x,y)[0]
    return res


def test_genepair(gene1, gene2, g1_new_alt, g1_new_ref, g2_new_alt, g2_new_ref, cell_express_filter=10):
    
    n = compute_correlation(g1_new_alt, g2_new_alt)
    n2 = compute_correlation(g1_new_ref, g2_new_ref)
    x1 = compute_correlation(g1_new_alt, g2_new_ref)
    x2 = compute_correlation(g1_new_ref, g2_new_alt)
    nonallele = compute_correlation(g1_new_ref+g1_new_alt, g2_new_ref+g2_new_alt)
    nmx = (n+n2)/2. - (x1+x2)/2.
    return gene1, gene2, n, n2, x1, x2, nmx, nonallele

if not o.load:
    
    new_alt = pd.read_hdf('alt_new.h5', 'alt_new')
    new_ref = pd.read_hdf('ref_new.h5', 'ref_new')
    new_altX = new_alt.values
    new_refX = new_ref.values
        
    
    # filter cells and genes.
    print('optional filtering to most highly epxressed genes, detected with new RNA in more than %i cells.' % cellt)
    genes_to_keep = []
    for gidx in range(len(new_alt.index.values)):
        tot = new_altX[gidx,:] + new_refX[gidx,:]
        if sum([t > 0 and 1 or 0 for t in tot]) >= cellt:
            genes_to_keep.append(new_alt.index[gidx])
    
    print('will keep %i genes.' % len(genes_to_keep))
    new_alt = new_alt.filter(items=genes_to_keep, axis=0)
    new_ref = new_ref.filter(items=genes_to_keep, axis=0)
    new_altX = new_alt.values
    new_refX = new_ref.values
    
    
    # read and index genes by starting position
    genes_d = {}
    genes_d_X = {}
    for line in open('genomic_positions.txt','r'):
        parts = line.strip().split('\t')
        if parts[0] == 'chr': continue
        gene = parts[0]
        chrom = parts[1]
        strand = parts[5]
        if strand == '1':
            start = int(parts[4])
        else:
            start = int(parts[2])
            
        if chrom == 'X':
            genes_d_X[gene] = [chrom, start, strand]
        else:
            genes_d[gene] = [chrom, start, strand]
    
    
    bins = [[],[],[],[],[]]
    separators = [100_000, 500_000, 1_500_000, 2_500_000, 3_500_000]
    
    genes = [gene.split('.')[0] for gene in new_alt.index.values]
    
    for gidx1, gene1 in enumerate(tqdm(genes)):
        for gidx2, gene2 in enumerate(genes):
            if gene1 not in genes_d or gene2 not in genes_d: continue
            if gidx1 <= gidx2: continue
            if genes_d[gene1][0] != genes_d[gene2][0]:
                continue
    
            # require high expression, at least 100 cells with expression
            
            if abs(genes_d[gene1][1]-genes_d[gene2][1]) < max(separators):
                sep_idx = 0
                dist = genes_d[gene1][1]-genes_d[gene2][1]
                while dist > separators[sep_idx]:
                    sep_idx += 1
                bins[sep_idx].append( (gene1, gene2, gidx1, gidx2) )
    
    binsX = [[],[],[],[],[]]
    for gidx1, gene1 in enumerate(tqdm(genes)):
        for gidx2, gene2 in enumerate(genes):
            if gene1 not in genes_d_X or gene2 not in genes_d_X: continue
            if gidx1 <= gidx2: continue
            if genes_d_X[gene1][0] != genes_d_X[gene2][0]:
                continue
    
            # require high expression, at least 100 cells with expression                                                                                                                                                                                                                           
    
            if abs(genes_d_X[gene1][1]-genes_d_X[gene2][1]) < max(separators):
                sep_idx = 0
                dist = genes_d_X[gene1][1]-genes_d_X[gene2][1]
                while dist > separators[sep_idx]:
                    sep_idx += 1
                binsX[sep_idx].append( (gene1, gene2, gidx1, gidx2) )
    
    
    print('gene pairs in bins')
    for s in (0,1,2,3,4):
        print(' < %i : %i   %i' % (separators[s], len(bins[s]), len(binsX[s])))
        print(bins[s][0])
                
    
    # now ready to read expression levels for genes on the chromosome
    nb_genes, nb_cells = new_alt.shape
    new_alt.set_axis(genes, axis=0)
    new_ref.set_axis(genes, axis=0)
    cells = new_alt.columns.values    
    
    results_nmean_X = []
    results_xmean_X = []
    results_nmx_X = []
    results_na_X = []
    
    results_nmean = []
    results_xmean = []
    results_nmx = []
    results_na = []
    
    for s in (0,1,2,3,4):
        res = Parallel(n_jobs=50)(delayed(test_genepair)(gene1, gene2, np.array(new_altX[g1_idx,:]), np.array(new_refX[g1_idx,:]),np.array(new_altX[g2_idx,:]), np.array(new_refX[g2_idx,:])) for gene1, gene2, g1_idx, g2_idx in bins[s])
    
        tmp_n = []
        tmp_x = []
        tmp_nmx = []
        tmp_na = []
        
        for vals in res:
            g1, g2, n, n2, x1, x2, nmx, na = vals
            
            if np.isnan(n): n=0
            if np.isnan(n2): n2=0
            if np.isnan(x1): x1=0
            if np.isnan(x2): x2=0
            if np.isnan(nmx): nmx=0
            if np.isnan(na): na = 0
            
            nmean = (n+n2)/2.
            xmean = (x1+x2)/2.
            if same_strand and genes_d[g1][2] != genes_d[g2][2]:
                continue
            
            tmp_n.append(nmean)
            tmp_x.append(xmean)
            tmp_nmx.append(nmx)
            tmp_na.append(na)
            
        results_nmean.append(tmp_n)
        results_xmean.append(tmp_x)
        results_nmx.append(tmp_nmx)
        results_na.append(tmp_na)
    
    for s in (0,1,2,3,4):
        res = Parallel(n_jobs=50)(delayed(test_genepair)(gene1, gene2, np.array(new_altX[g1_idx,:]), np.array(new_refX[g1_idx,:]),np.array(new_altX[g2_idx,:]), np.array(new_refX[g2_idx,:])) for gene1, gene2, g1_idx, g2_idx in binsX[s])
    
        tmp_n = []
        tmp_x = []
        tmp_nmx = []
        tmp_na = []
    
        for vals in res:
            g1, g2, n, n2, x1, x2, nmx, na = vals
    
            if np.isnan(n): n=0
            if np.isnan(n2): n2=0
            if np.isnan(x1): x1=0
            if np.isnan(x2): x2=0
            if np.isnan(nmx): nmx=0
            if np.isnan(na): na = 0
    
            nmean = (n+n2)/2.
            xmean = (x1+x2)/2.
            if same_strand and genes_d_X[g1][2] != genes_d_X[g2][2]:
                continue
    
            tmp_n.append(nmean)
            tmp_x.append(xmean)
            tmp_nmx.append(nmx)
            tmp_na.append(na)
    
        results_nmean_X.append(tmp_n)
        results_xmean_X.append(tmp_x)
        results_nmx_X.append(tmp_nmx)
        results_na_X.append(tmp_na)
    
    
    
    
    
    # https://stackoverflow.com/questions/2960864/how-to-save-all-the-variables-in-the-current-python-session
    import shelve
    filename='shelve.out'
    my_shelf = shelve.open(filename,'n') # 'n' for new
    for key in dir():
        try:
            my_shelf[key] = globals()[key]
        except TypeError:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
            print('ERROR1 shelving: {0}'.format(key))
        except AttributeError:
            print('ERROR2 shelving: {0}'.format(key))
    my_shelf.close()
else:
    import shelve
    filename='shelve.out'
    my_shelf = shelve.open(filename)
    for key in my_shelf:
        globals()[key]=my_shelf[key]
    my_shelf.close()

#  create composite figure

def melt_and_explode(coldict):
    df1 = pd.DataFrame.from_dict(coldict, orient='index').T
    df2 = df1.melt().dropna()
    df2['dist'] = np.tile(['<0.1', '<0.5', '<1.5', '<2.5', '<3.5'], len(df2)//5)
    print(df2)
    df3 = df2.explode('value').reset_index(drop=True)
    print(len(df3))
    df3['value'] = df3['value'].astype(float)
    return df3.dropna()

red_circle = dict(markerfacecolor='red', marker='.', markeredgecolor='none', markersize= 2)

sns_styleset()
plt.figure(figsize=(4,8))

plt.subplot(2,1,1)
autosomal_efg = melt_and_explode({'without_allelic':results_na, 'same_allele':results_nmean, 'other_allele':results_xmean})

'''
for var, ls in zip(autosomal_efg['variable'].unique(), ['-', '--', ':']):
    sns.kdeplot(x='value', data=autosomal_efg[autosomal_efg['variable']==var], hue='dist', linestyle=ls)
plt.xlabel('Spearman correlation')
'''
sns.boxplot(data=autosomal_efg, x='dist', y='value', hue='variable', flierprops=red_circle)
plt.ylabel('Spearman correlation difference (same - other)')
plt.xlabel('Genomic distance (Mbp)')
plt.gca().axhline(0)

plt.subplot(2,1,2)
difference_hg = melt_and_explode({'autosomal':results_nmx, 'Xchrom':results_nmx_X})
sns.boxplot(data=difference_hg, x='dist', y='value', hue='variable', flierprops=red_circle)
plt.ylabel('Spearman correlation difference (same - other)')
plt.xlabel('Genomic distance (Mbp)')
plt.gca().axhline(0)

sns.despine(offset={'left':3, 'bottom':0}, trim=True)
#plt.gcf().align_labels()
plt.tight_layout()


plt.savefig('figure2.png', dpi=150)
plt.savefig('figure2.pdf', dpi=300)

