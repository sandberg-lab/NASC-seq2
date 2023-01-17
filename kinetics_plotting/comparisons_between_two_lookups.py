import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import argparse, pandas, seaborn, numpy
from matplotlib import pyplot
from statsmodels.nonparametric import smoothers_lowess
from scipy import stats


def combinedP(pvalues, weights=None):
	""" takes arrays of p-values (preferably 1-sided) and weights (sample sizes if equal variance), returns p-value (1-sided if 1-sided input) """
	from math import sqrt
	Zs = [PtoZ(p) for p in pvalues]
	if weights is None:
		Zcombined = sum([Z for Z in Zs])/sqrt(len(Zs))
	else:
		Zcombined = sum([w*Z for w,Z in zip(weights, Zs)])/sqrt(sum([w**2 for w in weights]))	
	return ZtoP(Zcombined)
def PtoZ(p):
	from scipy.special import erfinv
	from math import sqrt
	return sqrt(2.0)*erfinv(2.0*p-1.0)
def ZtoP(Z):
	from scipy.special import erf
	from math import sqrt
	return 0.5*(1.0+erf(Z/sqrt(2.0)))

def load(path):
    if path.endswith('.csv'):
        df = pandas.read_csv(path, index_col=0)
        df['burstsize'] = df['burst_size']
        del df['burst_size']
        df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
        df['exprpct'] = df['expressionrate'].rank(pct=True, ascending=False)
        df['exprrank'] = df['expressionrate'].rank(pct=False, ascending=False)
        if not o.allML: df = df[df['success']]
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
        df['log10_%s'%(param)] = numpy.log10(df[param])
        
    if o.cut_at_index_dot:
        df.index = df.index.str.split('.').str[0]
    return df

def across_halves(df, o, params_here, res_here, descr):
    pyplot.clf()
    pyplot.figure(figsize=(8,8))
    drew_something = False
    for pi, xparam in enumerate(params_here):
        for ri, yparam in enumerate(res_here):
            if all(param in ('kon_x_degr', 'kon_/_degr') for param in (xparam, yparam)): continue
            pyplot.subplot(len(params_here), len(res_here), 1+ri+pi*len(res_here))
            xcol = (xparam if o.avoid_log else 'log10_'+xparam) + '_'+o.names[0]
            if xcol not in df:
                xcol = xparam
                if xcol not in df:
                    continue
            ycol = (yparam if o.avoid_log else 'log10_'+yparam) + '_'+o.names[1]
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
            
            if (xparam in ('kon', 'koff', 'ksyn') and yparam in ('kon', 'koff', 'ksyn')) or xparam==yparam:
                pyplot.gca().axline((0,0), (1,1))
            
    if drew_something:
        pyplot.suptitle('%i genes'%len(df))
        pyplot.tight_layout()
        pyplot.savefig('%s_%s_acrosshalves.pdf'%(o.prefix_out, descr))

def iter_dataframe(df):
    return ((idx, col, val) for col in df for idx, val in df[col].items())

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('table1')
    parser.add_argument('table2')
    parser.add_argument('prefix_out')
    parser.add_argument('--skip_borders', metavar=('kon_min', 'kon_max', 'koff_min', 'koff_max', 'ksyn_min', 'ksyn_max'), nargs=6, type=float)
    parser.add_argument('--expression_level_basis', nargs=2)
    parser.add_argument('--degradation', nargs=2)
    parser.add_argument('--filter_by_CI', nargs=6, metavar=('CIwidth', 'minratio', 'maxratio', 'maxfail%', 'num_filter_params', 'first_half_only'), help='e.g. 90 1 5 50 4 True or 50 1 2 50 2 False')
    parser.add_argument('--avoid_log', action='store_true')
    parser.add_argument('--combineP', action='store_true')
    parser.add_argument('--allML', action='store_true')
    parser.add_argument('--cut_at_index_dot', action='store_true')
    parser.add_argument('--genes_to_use')
    parser.add_argument('--genes_not_to_use')
    parser.add_argument('--names', nargs=2, default=['1', '2'])
    o = parser.parse_args()
    
    df1 = load(o.table1)
    df1.columns = [col+'_'+o.names[0] for col in df1.columns]
    df2 = load(o.table2)
    df2.columns = [col+'_'+o.names[1] for col in df2.columns]
    dfM = df1.merge(df2, left_index=True, right_index=True)
    print(len(df1), len(df2), len(dfM))
    if o.expression_level_basis is not None:
        for i, path in enumerate(o.expression_level_basis):
            dfE = load(path)
            dfE['expression'] = dfE['binary_avg']*dfE['nonzero_avg']
            dfE.columns = [col+'_'+o.names[i] for col in dfE.columns]
            for col in dfE.columns:
                if col in dfM.columns:
                    del dfE[col]
            dfM = dfM.merge(dfE, left_index=True, right_index=True)
            print(len(dfM))
    if o.degradation is not None:
        for i, path in enumerate(o.degradation):
            dfE = pandas.read_csv(path, index_col=0)
            dfE.columns = [col+'_'+o.names[i] for col in dfE.columns]
            for col in dfE.columns:
                if col in dfM.columns:
                    del dfE[col]
            dfM = dfM.merge(dfE, left_index=True, right_index=True)
            dfM['log10_kon_x_degr_'+o.names[i]] = numpy.log10(dfM['kon_'+o.names[i]] * dfM['degradationrate_'+o.names[i]])
            dfM['log10_kon_/_degr_'+o.names[i]] = numpy.log10(dfM['kon_'+o.names[i]] / dfM['degradationrate_'+o.names[i]])
        print(len(dfM))
    
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
    
    '''pyplot.clf()
    for param in ['kon', 'koff', 'ksyn', 'degradationrate']:
        xcol = param if o.avoid_log else 'log10_%s'%param
        if xcol not in dfM: continue
        seaborn.kdeplot(x=xcol, data=dfM)
    pyplot.tight_layout()
    pyplot.savefig(o.prefix_out + '_distr.pdf')'''
    
    if len(dfM)==0: raise Exception("No genes after filtering")
    
    across_halves(dfM, o, ['kon', 'koff', 'ksyn'], ['burstsize', 'expressionrate', 'switchingrate'], '1')
    across_halves(dfM, o, ['kon', 'koff', 'ksyn'], ['kon', 'koff', 'ksyn'], '2')
    across_halves(dfM, o, ['switchingrate', 'burstsize', 'expressionrate'], ['switchingrate', 'burstsize', 'expressionrate'], '3')
    across_halves(dfM, o, ['binary_avg', 'nonzero_avg', 'nonzero_cov'], ['kon', 'koff', 'ksyn'], '4')
    if o.names == ['new', 'total']: across_halves(dfM, o, ['kon_x_degr', 'kon', 'kon_/_degr'], ['kon_x_degr', 'kon', 'kon_/_degr'], '5')
    
    pyplot.clf()
    params = ['meanoccupancy', 'switchingrate', 'kon', 'binary_avg', 'expression', 'expressionrate', 'nonzero_avg', 'burstsize', 'ksyn', 'nonzero_cov', 'koff', 'degradationrate']
    if o.degradation is None or o.degradation[0] == o.degradation[1]: params.remove('degradationrate')
    params = [param for param in params if param+'_'+o.names[0] in dfM and param+'_'+o.names[1] in dfM]
    correlationmatrix = pandas.DataFrame(index=params, columns=params)
    pvaluematrix = pandas.DataFrame(index=range(len(params)), columns=range(len(params)))
    for xi, xparam in enumerate(params):
        for yi, yparam in enumerate(params):
            xcol = xparam+'_'+o.names[0]
            ycol = yparam+'_'+o.names[1]
            dfR = dfM[[xcol, ycol]].dropna()
            correlationmatrix.loc[yparam, xparam], pvaluematrix.loc[xi, yi] = stats.spearmanr(dfR[xcol], dfR[ycol])
    seaborn.heatmap(correlationmatrix.astype(float), cmap=seaborn.color_palette("vlag", as_cmap=True), vmin=-1, vmax=1, center=0, annot=True, square=True, xticklabels=True, yticklabels=True, fmt='.2f')
    for x,y,v in iter_dataframe(pvaluematrix):
        if x == y:
            pyplot.plot([x,x+1], [y,y], 'k-')
            pyplot.plot([x,x+1], [y+1,y+1], 'k-')
            pyplot.plot([x,x], [y,y+1], 'k-')
            pyplot.plot([x+1,x+1], [y,y+1], 'k-')
        elif o.combineP:
            if v > 0.05: print(v, pvaluematrix.loc[y, x], combinedP([v, pvaluematrix.loc[y, x]]))
            v = max([v, pvaluematrix.loc[y, x]])
        if v > 0.05:
            pyplot.text(x+0.5, y+0.75, 'ns', color='k', ha='center', va='center')
        
    pyplot.xlabel(o.names[0])
    pyplot.ylabel(o.names[1])
    pyplot.suptitle('%i genes'%len(dfM))
    pyplot.tight_layout()
    pyplot.savefig('%s_spearman.pdf'%o.prefix_out)
    