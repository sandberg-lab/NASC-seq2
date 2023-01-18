import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import sys
import numpy as np

adata = ad.read(sys.argv[1])
print(adata)

sums = adata.X.toarray().sum(axis=0) + 1
#print(adata.sum(axis=0))
plt.hist(np.log10(sums), bins=50)
plt.savefig('tmp.png')
plt.clf()

plt.subplot(211)
sumsN = adata.layers['REF_new'].toarray().sum(axis=0) + 1
#print(adata.sum(axis=0))                                                                                                                                                  
plt.hist(np.log10(sumsN), bins=50)
#plt.savefig('tnew.png')
plt.title('new REF')
plt.subplot(212)
sumsN = adata.layers['ALT_new'].toarray().sum(axis=0) + 1
#print(adata.sum(axis=0))
plt.hist(np.log10(sumsN), bins=50)
plt.title('new ALT')
plt.savefig('tnew.png')



plt.figure()
plt.scatter(np.log10(sums), np.log10(sumsN))
plt.savefig('scatter.png')


# cut-offs to use
# cells should have > 10k UMIs and >1k new UMIs

cells_to_use = []
for cell, tsum, nsum in zip(adata.var.index, sums, sumsN):
    if tsum > 10_000 and nsum > 1_000:
        cells_to_use.append(cell)

with open('cells_to_use.txt','w') as of:
    of.write("%s\n" % "\n".join(cells_to_use))

# genes to use, should have RNA in 1% of cells
genesums = (adata.X.toarray() > 0).sum(axis=1)
genes_to_use = []
for gene, tsum in zip(adata.obs.index, genesums):
    if tsum > 100:
        genes_to_use.append(gene)

with open('genes_to_use.txt','w') as of:
    of.write("%s\n" % "\n".join(genes_to_use))


# save data to h5py files
alt_new = adata.layers['ALT_new'].toarray()
df_alt_new = pd.DataFrame(alt_new, columns = adata.var.index, index= adata.obs.index)
print(df_alt_new.shape)
df_alt_new = df_alt_new[cells_to_use]
print(df_alt_new.shape)
df_alt_new = df_alt_new.filter(items=genes_to_use, axis=0)
print(df_alt_new.shape)
df_alt_new.to_hdf('alt_new.h5',key='alt_new', mode='w')


ref_new = adata.layers['REF_new'].toarray()
df_ref_new = pd.DataFrame(ref_new, columns = adata.var.index, index= adata.obs.index)
print(df_ref_new.shape)
df_ref_new = df_ref_new[cells_to_use]
print(df_ref_new.shape)
df_ref_new = df_ref_new.filter(items=genes_to_use, axis=0)
print(df_ref_new.shape)
df_ref_new.to_hdf('ref_new.h5',key='ref_new', mode='w')
df_ref_new.to_csv('ref_new.csv')
