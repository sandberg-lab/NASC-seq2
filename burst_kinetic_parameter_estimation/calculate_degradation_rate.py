import argparse, pandas, math, functools

def get_degradation_rate(time, fraction_new):
    # math.log = ln
    if fraction_new == 1:
        return float('inf')
    else:
        return -math.log(1-fraction_new)/time

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('newrna_csv')
    parser.add_argument('oldrna_csv')
    parser.add_argument('time', type=float)
    parser.add_argument('out_csv')
    parser.add_argument('--out_avg', default='/dev/stdout')
    parser.add_argument('--bootstrap', type=int, const=50, nargs='?')
    o = parser.parse_args()
    
    df_new = pandas.read_csv(o.newrna_csv, index_col=0)
    df_total = df_new + pandas.read_csv(o.oldrna_csv, index_col=0)
    
    ratios = (df_new.sum(axis=1) / df_total.sum(axis=1)).to_frame('fractionnew')
    
    with open(o.out_avg, 'wt') as outfh:
        fraction_new = float(df_new.sum().sum()/df_total.sum().sum())
        print('degradationrate:', get_degradation_rate(o.time, fraction_new), 'fractionnew:', fraction_new, file=outfh)
    
    ratios['degradationrate'] = ratios['fractionnew'].apply(functools.partial(get_degradation_rate, o.time))
    
    if o.bootstrap is not None:
        import random, numpy
        random.seed(10)
        boots = pandas.DataFrame()
        for bi in range(o.bootstrap):
            samples = random.choices(df_new.columns, k=len(df_new.columns))
            boots[bi] = (df_new[samples].sum(axis=1) / df_total[samples].sum(axis=1)).apply(functools.partial(get_degradation_rate, o.time))
        ratios[['degradationrate_95%conf_low', 'degradationrate_90%conf_low', 'degradationrate_50%conf_low', 'degradationrate_bootstrapmedian', 'degradationrate_50%conf_high', 'degradationrate_90%conf_high', 'degradationrate_95%conf_high']] = boots.apply(functools.partial(numpy.percentile, q=[2.5, 5, 25, 50, 75, 95, 97.5]), axis=1, result_type='expand')
    
    ratios.to_csv(o.out_csv)