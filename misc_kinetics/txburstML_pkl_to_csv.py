import pandas, argparse

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('pkl_in', metavar='pkl_from_txburstML')
    parser.add_argument('csv_out', nargs='?')
    o = parser.parse_args()
    
    if o.csv_out is None:
        o.csv_out = o.pkl_in.rsplit('.',1)[0] + '.csv'
    
    data = pandas.read_pickle(o.pkl_in)
    df_out = pandas.DataFrame(index=data.index)
    df_out.index.name = None
    df_out['kon'] = [V[0] for V in data[0]]
    df_out['koff'] = [V[1] for V in data[0]]
    df_out['ksyn'] = [V[2] for V in data[0]]
    df_out['burst_size'] = df_out['ksyn']/df_out['koff']
    df_out['success'] = data[1]
    
    df_out.to_csv(o.csv_out)
