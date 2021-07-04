#!/usr/bin/env python3


"""data
Usage: 03.deg.py [-d anndata] [--c1 condition1] [--c2 condition2] [-s subset_name] [-o obs_column_name] [-h help]


Options:
    -h help
    -d=anndatafile anndata obj with raw expression values
    --c1=condition 1 name ideally wildtype
    --c2=condition 2 name ideally mutant
    -s=cluster/group of cells to subset (eg. HSCs)
    -o=obs name that cluster falls under (eg. niki_label)
"""


import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import scipy as sp
import anndata
from statsmodels.stats import multitest as mt
import pandas as pd
sc.set_figure_params(dpi=90)
import pickle
import sklearn
import scipy.stats as st
from matplotlib.backends.backend_pdf import PdfPages
from docopt import docopt

from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list(name='gene_cmap',
                                        colors=[(0, 'lightgrey'), (0.1, 'yellow'), (0.5, 'red'), (1, 'darkred')])

anndatafile = '../ANALYSIS1/cluster_nikiLabelled.h5ad'
cond1 = 'WT'
cond2 = 'Tet2Hom'

arguments = docopt(__doc__)
#print(arguments)
anndatafile = arguments['-d']
cond1=arguments['--c1']
cond2=arguments['--c2']
clus = arguments['-s']
obs_name = arguments['-o']

def filter_genes(ad, x):
    boo = np.sum(ad.X > x, axis=0) > 0
    print('Total number of genes passing filtering is:', sum(boo))
    ad = ad[:, boo]
    return ad


def remove_zero_genes(d1, d2):
    temp = np.sum(d1.X, axis = 0)
    temp1 = np.sum(d2.X, axis = 0)
    cwt = pd.DataFrame({'genes': d1.var_names, 'wt': temp}, index=d1.var_names)
    chom = pd.DataFrame({'genes': d2.var_names, 'hom': temp1}, index = d2.var_names)
    counts = pd.concat([cwt, chom], axis = 1, join = 'inner' )
    cgenes = list(counts[(counts['wt'] > 10) & (counts['hom'] > 10)].index )
    return d1[:, cgenes], d2[:, cgenes]

def subset_to_clus(data, clus_list, name):
    boo = [clus in clus_list for clus in data.obs[name]]
    data= data[boo,:]
    return data

def remove_zero_genes(data1, data2):
    #remove genes which do nothing
    boo = [sum(data1.X[:,i]) == 0 and sum(data2.X[:,i]) == 0 for i in range(len(data1.var_names))]
    data1 = data1[:,np.invert(boo)]
    data2 = data2[:,np.invert(boo)]
    return data1, data2

def retain_shared_genes(data1, data2):
    boo1 = [name in data2.var_names for name in data1.var_names]
    boo2 = [name in data1.var_names for name in data2.var_names]
    data1 = data1[:,boo1]
    data2 = data2[:,boo2]
    print(len(data1.var_names))
    print(len(data2.var_names))
    print(data1.var_names == data2.var_names)
    return data1, data2


adata = sc.read(anndatafile)
print("Normalising and logging anndata")
sc.pp.normalize_total(adata, exclude_highly_expressed=True, target_sum=10000)

print(f"Subsetting to {clus}")
adata = subset_to_clus(adata, clus, obs_name)

exp_wt = adata[adata.obs['condition']=='WT',:].copy()
exp_hom = adata[adata.obs['condition']=='Tet2Hom',:].copy()
exp_wt.X = exp_wt.X.A
exp_hom.X = exp_hom.X.A


print(f"Shape of {cond1}")
print(exp_wt.shape)
print(f"Shape of {cond2}")
print(exp_hom.shape)


print("Filtering for decent genes")
# Filter genes
exp_wt = filter_genes(exp_wt, 3)
exp_hom = filter_genes(exp_hom, 3)

#subset genes in each experiment to only shared genes
exp_wt, exp_hom = retain_shared_genes(exp_wt, exp_hom)

print(f"Shape of {cond1}")
print(exp_wt.shape)
print(f"Shape of {cond2}")
print(exp_hom.shape)

print('Running scanpy DEG')
exp_adata = exp_wt.concatenate(exp_hom)
sc.pp.log1p(exp_adata)

sc.tl.rank_genes_groups(exp_adata, groupby='condition', reference = cond1, method='wilcoxon')
df = sc.get.rank_genes_groups_df(exp_adata, group = 'Tet2Hom')

print("Pickling for later")
bn = cond1 + "_" + cond2 + "_" + clus 
df.to_csv(bn + '.deg.csv')
df.to_pickle(bn + '.deg.pkl')










