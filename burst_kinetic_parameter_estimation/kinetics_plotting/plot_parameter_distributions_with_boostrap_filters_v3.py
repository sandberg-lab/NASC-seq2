import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import argparse, pandas, seaborn, math, numpy, itertools
from matplotlib import pyplot, ticker

logtext = 'log10'
logfunc = numpy.log10

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='+')
    parser.add_argument('-o', '--plot_out', required=True)
    parser.add_argument('-n', '--names', nargs='+')
    parser.add_argument('--degradation', nargs='+', default=[])
    parser.add_argument('--nrows', type=int)
    parser.add_argument('--drop_duplicates', choices=['no', 'yes', 'first', 'last'], default='no')
    parser.add_argument('--sep', default='\t')
    parser.add_argument('--line', nargs=2, action='append', metavar=('col', 'value'), default=[])
    parser.add_argument('-r', '--replratio_cutoff', type=float, nargs=3, default=[[float('-inf'), float('inf'), float('inf')]], action='append', metavar=('over', 'under', 'failmax%'))
    parser.add_argument('--conf', choices=['50', '90', '95'], default='95')
    parser.add_argument('--replmedians', action='store_true')
    parser.add_argument('--skip_borders', metavar=('kon_min', 'kon_max', 'koff_min', 'koff_max', 'ksyn_min', 'ksyn_max'), nargs=6, type=float)
    parser.add_argument('--conffilter', nargs='*', choices=['kon', 'koff', 'ksyn', 'burstsize'], default=[])
    parser.add_argument('--export_index')
    parser.add_argument('--single_plot', action='store_true')
    parser.add_argument('--same_genes', action='store_true')
    parser.add_argument('--expression_level_basis', nargs='+', default=[])
    parser.add_argument('--genes_to_use')
    parser.add_argument('--genes_not_to_use')
    parser.add_argument('--bw', type=float)
    parser.add_argument('--use_columns', nargs='+')
    parser.add_argument('--cut_at_index_dot', action='store_true')
    o = parser.parse_args()
    
    if o.names is None: o.names = list(range(len(o.infile)))
    
    dfs = []
    for path, name, expression_level_basis, degradation_file in itertools.zip_longest(o.infile, o.names, o.expression_level_basis, o.degradation):
        df = pandas.read_table(path, index_col=0, nrows=o.nrows, sep=o.sep)
        df = df.dropna(subset=['kon', 'koff', 'ksyn'])
        if o.cut_at_index_dot:
            df.index = df.index.str.split('.').str[0]
        df['expressionrate'] = df['ksyn'] * df['kon'] / (df['kon'] + df['koff'])
        df['meanoccupancy'] = df['kon'] / (df['kon'] + df['koff'])
        df['switchingrate'] = 1 / (1/df['kon'] + 1/df['koff'])
        if o.drop_duplicates != 'no':
            df = df.drop_duplicates(subset=['kon', 'koff', 'ksyn'], keep=(False if o.drop_duplicates=='yes' else o.drop_duplicates))
        if o.skip_borders is not None:
            for pi, param in enumerate(['kon', 'koff', 'ksyn']):
                minv, maxv = o.skip_borders[2*pi: 2*pi+2]
                col = param+'_bootstrapmedian' if o.replmedians else param
                df = df[(df[col] > minv) & (df[col] < maxv)]
        if o.genes_to_use is not None:
            with open(o.genes_to_use, 'rt') as infh:
                genes = set(line.strip() for line in infh)
            df = df.loc[df.index.str.split('.').str[0].isin(genes) | df.index.isin(genes), :]
        if o.genes_not_to_use is not None:
            with open(o.genes_not_to_use, 'rt') as infh:
                genes = set(line.strip() for line in infh)
            df = df.loc[~(df.index.str.split('.').str[0].isin(genes) | df.index.isin(genes)), :]
        if expression_level_basis is not None:
            dfE = pandas.read_table(expression_level_basis, index_col=0)
            if o.cut_at_index_dot:
                dfE.index = df.index.str.split('.').str[0]
            for col in dfE.columns:
                if col in df.columns:
                    del dfE[col]
            df = df.merge(dfE, left_index=True, right_index=True)
        if degradation_file is not None:
            dfE = pandas.read_csv(degradation_file, index_col=0)
            if o.cut_at_index_dot:
                dfE.index = df.index.str.split('.').str[0]
            for col in dfE.columns:
                if col in df.columns:
                    del dfE[col]
            df = df.merge(dfE, left_index=True, right_index=True)
        df.index.name = name
        for col in o.conffilter:
            repl_ratio = df[col+'_'+o.conf+'%conf_high'] / df[col+'_'+o.conf+'%conf_low']
            lim_low, lim_high, max_fail = o.replratio_cutoff[-1]
            if col == 'degradationrate':
                df = df[(repl_ratio > lim_low) & (repl_ratio < lim_high)]
            else:
                df = df[(repl_ratio > lim_low) & (repl_ratio < lim_high) & (df[col+'_bootstrapfail%'] < max_fail)]
        dfs.append(df)
        #df.columns = [col+'_'+name for col in df.columns]
        #dfM = df if dfM is None else dfM.merge(df, left_index=True, right_index=True)
    
    
    
    
    if o.same_genes:
        genelist = None
        columns = ['kon', 'koff', 'ksyn', 'burstsize'] if o.use_columns is None else o.use_columns
        for (df, lim, col), color in zip(itertools.product(dfs, o.replratio_cutoff[1:], columns), ['red', 'green', 'blue', 'grey', 'cyan', 'magenta', 'brown', 'purple', 'orange', 'black']):
            (lim_low, lim_high, max_fail) = lim
            df[logtext+'_'+col] = logfunc(df[col+'_bootstrapmedian']) if o.replmedians else logfunc(df[col])
            repl_ratio = df[col+'_'+o.conf+'%conf_high'] / df[col+'_'+o.conf+'%conf_low']
            if col == 'degradationrate':
                dfR = df[(repl_ratio > lim_low) & (repl_ratio < lim_high)]
            else:
                dfR = df[(repl_ratio > lim_low) & (repl_ratio < lim_high) & (df[col+'_bootstrapfail%'] < max_fail)]
            genelist = set(dfR.index) if genelist is None else genelist & set(dfR.index)
        dfs = [df.loc[genelist, :] for df in dfs]
    
    if o.export_index is not None:
        with open(o.export_index, 'wt') as outfh:
            for gene in dfs[0].index:
                print(gene, file=outfh)
                
    
    if o.single_plot:
        columns = o.use_columns if o.use_columns is not None else ['kon', 'koff', 'ksyn', 'burstsize', 'degradationrate'] if len(dfs)==1 else ['kon', 'ksyn', 'koff']
        for (df, linestyle) in zip(dfs, itertools.cycle(['-', '--', '-.', ':'])):
            for (lim, col), color in zip(itertools.product(o.replratio_cutoff[1:], columns), ['red', 'green', 'blue', 'grey', 'cyan', 'magenta', 'brown', 'purple', 'orange', 'black']):
                (lim_low, lim_high, max_fail) = lim
                df[logtext+'_'+col] = logfunc(df[col+'_bootstrapmedian']) if o.replmedians else logfunc(df[col])
                repl_ratio = df[col+'_'+o.conf+'%conf_high'] / df[col+'_'+o.conf+'%conf_low']
                if col == 'degradationrate':
                    dfR = df[(repl_ratio > lim_low) & (repl_ratio < lim_high)]
                else:
                    dfR = df[(repl_ratio > lim_low) & (repl_ratio < lim_high) & (df[col+'_bootstrapfail%'] < max_fail)]
                seaborn.kdeplot(x=col, log_scale=True, data=dfR, linestyle=linestyle, color=color, label='%s%s%s%s, %i genes'%(df.index.name+' ' if len(dfs)>1 else '',  col, ' (× yield)' if col in ('ksyn','burst_size','burstsize') else '', '/hour' if col in ('kon','koff','ksyn') else '', len(dfR)), bw_method=o.bw)
        for param, val in o.line:
            if param == col:
                pyplot.gca().axvline(logfunc(float(val)))
        from matplotlib import ticker
        pyplot.xlabel('')
        pyplot.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator())
        pyplot.legend()
    else:
        columns = ['kon', 'koff', 'ksyn', 'burstsize'] if o.use_columns is None else o.use_columns
        subplot_rows = int(math.floor(len(columns)**0.5))
        subplot_cols = int(math.ceil(len(columns)/subplot_rows))
        for col_i, col in enumerate(columns):
            pyplot.subplot(subplot_rows, subplot_cols, col_i+1)
            for (df, lim), color in zip(itertools.product(dfs, o.replratio_cutoff), ['red', 'green', 'blue', 'grey', 'cyan', 'magenta', 'brown', 'purple', 'orange', 'black']):
                (lim_low, lim_high, max_fail) = lim
                df[logtext+'_'+col] = logfunc(df[col+'_bootstrapmedian']) if o.replmedians else logfunc(df[col])
                repl_ratio = df[col+'_'+o.conf+'%conf_high'] / df[col+'_'+o.conf+'%conf_low']
                dfR = df[(repl_ratio > lim_low) & (repl_ratio < lim_high) & (df[col+'_bootstrapfail%'] < max_fail)]
                print(col, lim, len(dfR))
                seaborn.kdeplot(x=logtext+'_'+col, data=dfR, color=color, bw_method=o.bw)
            for param, val in o.line:
                if param == col:
                    pyplot.gca().axvline(logfunc(float(val)))
            pyplot.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
            if not o.conffilter: pyplot.ylim(0,1)
    pyplot.tight_layout()
    pyplot.savefig(o.plot_out)

