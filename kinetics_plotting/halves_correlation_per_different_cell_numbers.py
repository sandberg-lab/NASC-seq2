import argparse, pandas, seaborn, numpy, collections
from matplotlib import pyplot
from statsmodels.nonparametric import smoothers_lowess
from scipy import stats

def skip_borders(df):
    if o.skip_borders is not None:
        #df = df[~(df['kon'].isin(o.skip_borders[:2]) | df['koff'].isin(o.skip_borders[2:4]) | df['ksyn'].isin(o.skip_borders[4:6]))] # this was buggy before 3 oct 2022
        
        if o.border_proximity is None:
            try:
                df = df[(df['kon'] > o.skip_borders[0]) & (df['kon'] < o.skip_borders[1]) & (df['koff'] > o.skip_borders[2]) & (df['koff'] < o.skip_borders[3]) & (df['ksyn'] > o.skip_borders[4]) & (df['ksyn'] < o.skip_borders[5])]
            except:
                print(o.skip_borders)
                raise
        else:
            for param, border in zip(['kon', 'kon', 'koff', 'koff', 'ksyn', 'ksyn'], o.skip_borders):
                df = df[(df[param] > border*o.border_proximity) | (df[param] < border/o.border_proximity)]
    return df

def load(path, o, colsuffix):
    if path.endswith('.csv'):
        df = pandas.read_csv(path, index_col=0)
        df['burstsize'] = df['burst_size']
        del df['burst_size']
        df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
        df['exprpct'] = df['expressionrate'].rank(pct=True, ascending=False)
        df['exprrank'] = df['expressionrate'].rank(pct=False, ascending=False)
        df = df[df['success']]
    else:
        df = pandas.read_table(path, index_col=0)
        df['exprpct'] = (df['binary_avg'] * df['nonzero_avg']).rank(pct=True, ascending=False)
        df['exprrank'] = (df['binary_avg'] * df['nonzero_avg']).rank(pct=False, ascending=False)
        df = df.dropna(subset=['kon', 'koff', 'ksyn'] if o.dropna is None else o.dropna)
        
        if o.filter_by_CI is not None:
            lim_low, lim_high, max_fail = map(float, o.filter_by_CI[1:4])
            conf = o.filter_by_CI[0]
            for col in ['kon', 'burstsize', 'koff', 'ksyn'][:int(o.filter_by_CI[4])]:
                repl_ratio = df[col+'_'+conf+'%conf_high'] / df[col+'_'+conf+'%conf_low']
                df = df[(repl_ratio > lim_low) & (repl_ratio < lim_high) & (df[col+'_bootstrapfail%'] < max_fail)]
        
    df = skip_borders(df)
    
    df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
    df['meanoccupancy'] = df['kon'] / (df['kon'] + df['koff'])
    df['switchingrate'] = 1 / (1/df['kon'] + 1/df['koff'])
    
    if o.min_perc_cells is not None:
        df = df[df['binary_avg'] >= o.min_perc_cells/100]
    if o.min_kon is not None:
        df = df[df['kon'] >= o.min_kon]
    if o.min_bs is not None:
        df = df[df['burstsize'] >= o.min_bs]
    if o.max_koff is not None:
        df = df[df['koff'] < o.max_koff]
    if o.max_ksyn is not None:
        df = df[df['ksyn'] < o.max_ksyn]
    if o.genes_to_use is not None:
        with open(o.genes_to_use, 'rt') as infh:
            genes = set(line.strip() for line in infh)
        df = df.loc[df.index.str.split('.').str[0].isin(genes) | df.index.isin(genes), :]
    if o.genes_not_to_use is not None:
        with open(o.genes_not_to_use, 'rt') as infh:
            genes = set(line.strip() for line in infh)
        df = df.loc[~(df.index.str.split('.').str[0].isin(genes) | df.index.isin(genes)), :]
    
    
    df.columns = ['%s_%s'%(col, colsuffix) for col in df.columns]
    for param in ('kon', 'koff', 'ksyn', 'burstsize', 'expressionrate', 'meanoccupancy', 'switchingrate'):
        df['log2_%s_%s'%(param, colsuffix)] = numpy.log2(df['%s_%s'%(param, colsuffix)])
    if not path.endswith('.csv'):
        for param in ('kon', 'koff', 'ksyn', 'burstsize'):
            df['log2_CI_ratio_%s_%s'%(param, colsuffix)] = numpy.log2(df['%s_95%%conf_high_%s'%(param, colsuffix)]) - numpy.log2(df['%s_95%%conf_low_%s'%(param, colsuffix)])
    return df

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--files', nargs=3, action='append', metavar=('numcells', 'half1', 'half2'), required=True)
    parser.add_argument('-o', '--prefix_out', required=True)
    parser.add_argument('--skip_borders', metavar=('kon_min', 'kon_max', 'koff_min', 'koff_max', 'ksyn_min', 'ksyn_max'), nargs=6, type=float)
    parser.add_argument('--border_proximity', nargs='?', const=1.15, type=float)
    parser.add_argument('--min_perc_cells', type=float)
    parser.add_argument('--min_kon', type=float)
    parser.add_argument('--min_bs', type=float)
    parser.add_argument('--max_koff', type=float)
    parser.add_argument('--max_ksyn', type=float)
    parser.add_argument('--genes_to_use')
    parser.add_argument('--genes_not_to_use')
    parser.add_argument('--names', nargs=2, default=['1', '2'])
    parser.add_argument('--top_n_genes', nargs='+', type=int, default=[11313, 8000, 5656, 4000, 2828, 2000, 1414, 1000, 707, 500, 353])
    parser.add_argument('--filter_by_CI', nargs=5, metavar=('CIwidth', 'minratio', 'maxratio', 'maxfail%', 'num_filter_params'), help='e.g. 90 1 5 50 4 or 50 1 2 50 2')
    parser.add_argument('--dropna', nargs='*')
    o = parser.parse_args()
    
    cellnum_to_df = dict()
    for cellnum_str, file1, file2 in o.files:
        df1 = load(file1, o, o.names[0])
        df2 = load(file2, o, o.names[1])
        dfM = df1.merge(df2, left_index=True, right_index=True)
        cellnum_to_df[cellnum_str] = dfM
    '''
    pyplot.clf()
    for pi, param in enumerate(['koff', 'ksyn', 'kon', 'burstsize']):
        pyplot.subplot(2, 2, pi+1)
        heatmap_df = pandas.DataFrame()
        for exprpct_cutoff in (1, 0.7, 0.5, 0.35, 0.25, 0.18, 0.13, 0.09, 0.06):
            for cellnum_str, dfM in cellnum_to_df.items():
                dfR = dfM[dfM['exprpct_%s'%o.names[0]] <= exprpct_cutoff]
                V1 = dfR['log2_%s_%s'%(param, o.names[0])]
                V2 = dfR['log2_%s_%s'%(param, o.names[1])]
                r = ((V1-V2)/(V1+V2)).abs().median()
                heatmap_df.loc[cellnum_str, exprpct_cutoff] = r
        seaborn.heatmap(heatmap_df, vmin=0, center=0, cbar_kws={'label':'|v1-v2|/(v1+v2)_'+param}, annot=True, annot_kws={'size':4})
        pyplot.xlabel('Fraction of genes included')
        pyplot.ylabel('Number of cells per half')
    pyplot.tight_layout()
    pyplot.savefig(o.prefix_out + 'halveserror.pdf')
    '''
    
    '''
    pyplot.clf()
    seaborn.set(font_scale=0.7)
    for pi, param in enumerate(['koff', 'ksyn', 'kon', 'burstsize']):
        pyplot.subplot(2, 2, pi+1)
        heatmap_df = pandas.DataFrame()
        num_genes_per_cutoff = collections.defaultdict(list)
        for exprpct_cutoff in (1, 0.7, 0.5, 0.35, 0.25, 0.18, 0.13, 0.09, 0.06):
            for cellnum_str, dfM in cellnum_to_df.items():
                dfR = dfM[dfM['exprpct_%s'%o.names[0]] <= exprpct_cutoff]
                num_genes_per_cutoff[exprpct_cutoff].append(len(dfR))
                V1 = dfR['log2_%s_%s'%(param, o.names[0])]
                V2 = dfR['log2_%s_%s'%(param, o.names[1])]
                r, P = stats.pearsonr(V1, V2)
                
                heatmap_df.loc[cellnum_str, exprpct_cutoff] = r
        heatmap_df.columns = ['%.0f%% (avg %.0f)'%(100*exprpct_cutoff, numpy.mean(num_genes_per_cutoff[exprpct_cutoff])) for exprpct_cutoff in heatmap_df.columns]
        seaborn.heatmap(heatmap_df, vmin=-1, vmax=1, center=0, cbar_kws={'label':'r_halves_'+param}, annot=True, annot_kws={'size':4})
        pyplot.xlabel('Fraction of genes included')
        pyplot.ylabel('Number of cells per half')
        pyplot.gca().set_yticklabels(pyplot.gca().get_yticklabels(),rotation = 0)
    pyplot.tight_layout()
    pyplot.savefig(o.prefix_out + 'halvescorrelation.pdf')
    
    '''
    
    first_df = cellnum_to_df[o.files[0][0]]
    genes_per_exprrank_cutoff = {exprrank_cutoff : first_df.index[first_df['exprrank_%s'%o.names[0]] <= exprrank_cutoff] for exprrank_cutoff in o.top_n_genes}
    
    pyplot.clf()
    seaborn.set(font_scale=0.7)
    for pi, param in enumerate(['koff', 'ksyn', 'kon', 'burstsize']):
        pyplot.subplot(2, 2, pi+1)
        heatmap_df = pandas.DataFrame()
        pct_genes_per_cutoff = collections.defaultdict(list)
        #for exprrank_cutoff in (11313, 8000, 5656, 4000, 2828, 2000, 1414, 1000, 707, 500, 353):
        for exprrank_cutoff, genes in genes_per_exprrank_cutoff.items():
            for cellnum_str, dfM in cellnum_to_df.items():
                if cellnum_str == 'skip': continue
                #dfR = dfM[dfM['exprrank_%s'%o.names[0]] <= exprrank_cutoff]
                dfR = dfM.loc[dfM.index.intersection(genes), :]
                pct_genes_per_cutoff[exprrank_cutoff].append(dfR['exprpct_'+o.names[0]].max())
                V1 = dfR['log2_%s_%s'%(param, o.names[0])]
                V2 = dfR['log2_%s_%s'%(param, o.names[1])]
                try: r, P = stats.pearsonr(V1, V2)
                except ValueError: r = None
                heatmap_df.loc[cellnum_str, exprrank_cutoff] = r
        heatmap_df.columns = ['%i (avg %.0f%%)'%(exprrank_cutoff, 100*numpy.mean(pct_genes_per_cutoff[exprrank_cutoff])) for exprrank_cutoff in heatmap_df.columns]
        seaborn.heatmap(heatmap_df, vmin=-1, vmax=1, center=0, cbar_kws={'label':'r_halves_'+param}, annot=True, annot_kws={'size':4})
        pyplot.xlabel('Number of genes included')
        if not any('lookup' in row or 'ML' in row for row in heatmap_df.index):
            pyplot.ylabel('Number of cells per half')
        pyplot.gca().set_yticklabels(pyplot.gca().get_yticklabels(),rotation = 0)
    pyplot.tight_layout()
    pyplot.savefig(o.prefix_out + 'halvescorr_topngenes.pdf')