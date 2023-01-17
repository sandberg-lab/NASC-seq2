import argparse, pandas, seaborn, numpy
from matplotlib import pyplot
from statsmodels.nonparametric import smoothers_lowess
from scipy import stats

def load(path):
    if path.endswith('.csv'):
        df = pandas.read_csv(path, index_col=0)
        df['burstsize'] = df['burst_size']
        del df['burst_size']
        df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
        df['exprpct'] = df['expressionrate'].rank(pct=True, ascending=False)
        df['exprrank'] = df['expressionrate'].rank(pct=False, ascending=False)
        df = df[df['success']]
        if o.skip_borders is not None:
            df = df[(df['kon'] > o.skip_borders[0]) & (df['kon'] < o.skip_borders[1]) & (df['koff'] > o.skip_borders[2]) & (df['koff'] < o.skip_borders[3]) & (df['ksyn'] > o.skip_borders[4]) & (df['ksyn'] < o.skip_borders[5])]
    else:
        df = pandas.read_table(path, index_col=0)
        if o.skip_borders is not None:
            try:
                df = df[(df['kon'] > o.skip_borders[0]) & (df['kon'] < o.skip_borders[1]) & (df['koff'] > o.skip_borders[2]) & (df['koff'] < o.skip_borders[3]) & (df['ksyn'] > o.skip_borders[4]) & (df['ksyn'] < o.skip_borders[5])]
            except:
                print(o.skip_borders)
                raise
        filtercolumns = ['kon', 'burstsize', 'koff', 'ksyn']
        if o.filter_by_CI is not None:
            lim_low, lim_high, max_fail = map(float, o.filter_by_CI[1:4])
            conf = o.filter_by_CI[0]
            filtercolumns = ['kon', 'burstsize', 'koff', 'ksyn'][:int(o.filter_by_CI[4])] if int(o.filter_by_CI[4])>=0 else [['kon', 'burstsize', 'koff', 'ksyn'][-int(o.filter_by_CI[4])]]
            print(filtercolumns)
            for col in filtercolumns:
                repl_ratio = df[col+'_'+conf+'%conf_high'] / df[col+'_'+conf+'%conf_low']
                if col == 'degradationrate':
                    df = df[(repl_ratio > lim_low) & (repl_ratio < lim_high)]
                else:
                    df = df[(repl_ratio > lim_low) & (repl_ratio < lim_high) & (df[col+'_bootstrapfail%'] < max_fail)]
                print(len(df))
        df['exprpct'] = (df['binary_avg'] * df['nonzero_avg']).rank(pct=True, ascending=False)
        df['exprrank'] = (df['binary_avg'] * df['nonzero_avg']).rank(pct=False, ascending=False)
        df = df.dropna(subset=filtercolumns)
        print(len(df))
    
    
    df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
    df['meanoccupancy'] = df['kon'] / (df['kon'] + df['koff'])
    df['switchingrate'] = 1 / (1/df['kon'] + 1/df['koff'])
    for param in ('kon', 'koff', 'ksyn', 'burstsize', 'expressionrate', 'meanoccupancy', 'switchingrate'):
        df['log2_%s'%(param)] = numpy.log2(df[param])
    return df

def across_halves(df, o, params_here, res_here, descr):
    pyplot.clf()
    pyplot.figure(figsize=(8,8))
    drew_something = False
    for pi, xparam in enumerate(params_here):
        for ri, yparam in enumerate(res_here):
            if xparam == yparam: continue
            pyplot.subplot(len(params_here), len(res_here), 1+ri+pi*len(res_here))
            xcol = xparam if o.avoid_log else 'log2_'+xparam
            if xcol not in df:
                xcol = xparam
                if xcol not in df:
                    continue
            ycol = yparam if o.avoid_log else 'log2_'+yparam
            if ycol not in df:
                ycol = yparam
                if ycol not in df:
                    continue
            dfR = df[[xcol, ycol]].dropna()
            if len(dfR) == 0:
                continue
            drew_something = True
            r, P = stats.spearmanr(dfR[xcol], dfR[ycol])
            pyplot.title('Spearman r=%.2f'%r)
            #seaborn.scatterplot(x=xcol, y=ycol, data=dfR, s=5, linewidth=0, color='k', alpha=0.3)
            seaborn.kdeplot(x=xcol, y=ycol, data=dfR, fill=True)
            xycurve = smoothers_lowess.lowess(df[ycol], df[xcol])
            pyplot.plot(xycurve[:, 0], xycurve[:, 1], 'r--')
            
            if xcol.startswith('log') and ycol.startswith('log'):
                xlims = pyplot.gca().get_xlim()
                ylims = pyplot.gca().get_ylim()
                if xlims[1]-xlims[0] > ylims[1] - ylims[0]:
                    mid = (ylims[1] + ylims[0])/2
                    halfwidth = (xlims[1]-xlims[0])/2
                    pyplot.ylim(mid - halfwidth, mid + halfwidth)
                else:
                    mid = (xlims[1] + xlims[0])/2
                    halfwidth = (ylims[1]-ylims[0])/2
                    pyplot.xlim(mid - halfwidth, mid + halfwidth)
                pyplot.gca().set_aspect('equal')
            
            if xparam in ('kon', 'koff', 'ksyn') and yparam in ('kon', 'koff', 'ksyn'):
                pyplot.gca().axline((0,0), (1,1))
            
    if drew_something:
        pyplot.suptitle('%i genes'%len(df))
        pyplot.tight_layout()
        pyplot.savefig('%s_%s_acrosshalves.pdf'%(o.prefix_out, descr))

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('table')
    parser.add_argument('prefix_out')
    parser.add_argument('--skip_borders', metavar=('kon_min', 'kon_max', 'koff_min', 'koff_max', 'ksyn_min', 'ksyn_max'), nargs=6, type=float)
    parser.add_argument('--expression_level_basis')
    parser.add_argument('--filter_by_CI', nargs=6, metavar=('CIwidth', 'minratio', 'maxratio', 'maxfail%', 'num_filter_params', 'first_half_only'), help='e.g. 90 1 5 50 4 True or 50 1 2 50 2 False')
    parser.add_argument('--avoid_log', action='store_true')
    parser.add_argument('--genes_to_use')
    parser.add_argument('--genes_not_to_use')
    o = parser.parse_args()
    
    dfM = load(o.table)
    
    if o.expression_level_basis is not None:
        dfE = load(o.expression_level_basis)
        dfE['Index_of_dispersion'] = dfE['nonzero_cov']**2 * dfE['nonzero_avg'] # aka Fano factor
        for col in dfE.columns:
            if col in dfM.columns:
                del dfE[col] # was del dfM[col] before 20 oct 2022 so was bug
        dfM = dfM.merge(dfE, left_index=True, right_index=True)
    
    if o.genes_to_use is not None:
        with open(o.genes_to_use, 'rt') as infh:
            genes = set(line.strip() for line in infh)
        dfM = dfM.loc[dfM.index.str.split('.').str[0].isin(genes) | dfM.index.isin(genes), :]
        print(len(dfM))
    if o.genes_not_to_use is not None:
        with open(o.genes_not_to_use, 'rt') as infh:
            genes = set(line.strip() for line in infh)
        dfM = dfM.loc[~(dfM.index.str.split('.').str[0].isin(genes) | dfM.index.isin(genes)), :]
        print(len(dfM))
    
    pyplot.clf()
    for pi, (param, xlim) in enumerate((('kon', (-8, 2)), ('koff', (-1, 9)), ('ksyn', (0, 10)), ('burstsize', (-3, 7)))):
        pyplot.subplot(2,2,pi+1)
        xcol = param if o.avoid_log else 'log2_%s'%param
        seaborn.kdeplot(x=xcol, data=dfM, color='black')
        if not o.avoid_log:
            pyplot.xticks(list(range(xlim[0], xlim[1]+1)), ['%.2f'%(2**x) if x < 0 else str(2**x) for x in range(xlim[0], xlim[1]+1)])
            #pyplot.xlim(*xlim)
        pyplot.xlabel(param)
    pyplot.tight_layout()
    pyplot.savefig(o.prefix_out + '_distr.pdf')
    
    across_halves(dfM, o, ['kon', 'koff', 'ksyn'], ['burstsize', 'expressionrate', 'switchingrate'], '1')
    across_halves(dfM, o, ['kon', 'koff', 'ksyn'], ['kon', 'koff', 'ksyn'], '2')
    across_halves(dfM, o, ['switchingrate', 'burstsize', 'expressionrate'], ['switchingrate', 'burstsize', 'expressionrate'], '3')
    across_halves(dfM, o, ['binary_avg', 'nonzero_avg', 'nonzero_cov', 'Index_of_dispersion'], ['kon', 'koff', 'ksyn'], '4')