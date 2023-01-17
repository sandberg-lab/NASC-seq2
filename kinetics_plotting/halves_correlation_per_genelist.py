import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import argparse, pandas, seaborn, numpy, collections
from matplotlib import pyplot
from scipy import stats

def load(path, o, colsuffix):
    if path.endswith('.csv'):
        df = pandas.read_csv(path, index_col=0)
        df['burstsize'] = df['burst_size']
        del df['burst_size']
        df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
        df['exprpct'] = df['expressionrate'].rank(pct=True, ascending=False)
        df['exprrank'] = df['expressionrate'].rank(pct=False, ascending=False)
        #df = df[df['success']]
    else:
        df = pandas.read_table(path, index_col=0)
        df['exprpct'] = (df['binary_avg'] * df['nonzero_avg']).rank(pct=True, ascending=False)
        df['exprrank'] = (df['binary_avg'] * df['nonzero_avg']).rank(pct=False, ascending=False)
        #df = df.dropna(subset=['kon', 'koff', 'ksyn'])
        
    
    df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
    df['meanoccupancy'] = df['kon'] / (df['kon'] + df['koff'])
    df['switchingrate'] = 1 / (1/df['kon'] + 1/df['koff'])
    
    
    df.columns = ['%s_%s'%(col, colsuffix) for col in df.columns]
    for param in ('kon', 'koff', 'ksyn', 'burstsize', 'expressionrate', 'meanoccupancy', 'switchingrate'):
        df['log2_%s_%s'%(param, colsuffix)] = numpy.log2(df['%s_%s'%(param, colsuffix)])
    return df

def load_gene_list(path):
    with open(path, 'rt') as infh:
        return [line.strip() for line in infh]

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--files', nargs=3, action='append', metavar=('numcells', 'half1', 'half2'), required=True)
    parser.add_argument('-o', '--prefix_out', required=True)
    parser.add_argument('--names', nargs=2, default=['1', '2'])
    parser.add_argument('-g', '--genelist', nargs=2, required=True, action='append', metavar=('listname', 'path'))
    parser.add_argument('--palette', choices=['Spectral', 'vlag'], default='vlag')
    o = parser.parse_args()
    
    cellnum_to_df = dict()
    for cellnum_str, file1, file2 in o.files:
        df1 = load(file1, o, o.names[0])
        df2 = load(file2, o, o.names[1])
        dfM = df1.merge(df2, left_index=True, right_index=True)
        cellnum_to_df[cellnum_str] = dfM
    first_df = cellnum_to_df[o.files[0][0]]
    genes_per_exprrank_cutoff = {listname : load_gene_list(path) for listname, path in o.genelist}
    
    pyplot.clf()
    seaborn.set(font_scale=0.7)
    for pi, param in enumerate(['koff', 'ksyn', 'kon', 'burstsize']):
        pyplot.subplot(2, 2, pi+1)
        heatmap_df = pandas.DataFrame()
        num_genes_per_cutoff = collections.defaultdict(list)
        for exprrank_cutoff, genes in genes_per_exprrank_cutoff.items():
            for cellnum_str, dfM in cellnum_to_df.items():
                if cellnum_str == 'skip': continue
                #dfR = dfM[dfM['exprrank_%s'%o.names[0]] <= exprrank_cutoff]
                dfR = dfM.loc[dfM.index.str.split('.').str[0].isin(genes), :]
                num_genes_per_cutoff[exprrank_cutoff].append(len(dfR))
                V1 = dfR['%s_%s'%(param, o.names[0])]
                V2 = dfR['%s_%s'%(param, o.names[1])]
                try: r, P = stats.spearmanr(V1, V2)
                except ValueError: print('meh'); r = None
                heatmap_df.loc[cellnum_str, exprrank_cutoff] = r
        heatmap_df.columns = ['%s (%i genes)'%(col, int(numpy.mean(num_genes_per_cutoff[col]))) for col in heatmap_df.columns]
        seaborn.heatmap(heatmap_df, vmin=-1, vmax=1, center=0, cbar_kws={'label':'r_halves_'+param}, annot=True, annot_kws={'size':4}, cmap=seaborn.color_palette(o.palette, as_cmap=True), square=len(heatmap_df.index)==len(heatmap_df.columns))
        pyplot.gca().set_yticklabels(pyplot.gca().get_yticklabels(),rotation = 0)
    pyplot.tight_layout()
    pyplot.savefig(o.prefix_out + 'halvescorr_genelists.pdf')