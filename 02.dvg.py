#!/usr/bin/env python3


"""data
Usage: 02.dvg.py [--a1 anndata1] [--a2 anndata2] [--c1 condition1] [--c2 condition2] [-s subset_name] [-o obs_column_name] [-h help]


Options:
    -h help
    --a1=anndata1 (Quantile normalised expression values)
    --a2=anndata2 (Quantile normalised expression values)
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
file1=arguments['--a1']
file2=arguments['--a2']
cond1=arguments['--c1']
cond2=arguments['--c2']
clus = arguments['-s']
obs_name = arguments['-o']



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


# Reading in adata
data_wt = sc.read(file1)
data_hom = sc.read(file2)

# Subsetting
data_wt_st = subset_to_clus(data_wt, clus, obs_name)
data_hom_st = subset_to_clus(data_hom, clus, obs_name)


# Remove genes that do not have any data in them
print("Removing genes with no counts in them")
data_wt_st, data_hom_st = remove_zero_genes(data_wt_st, data_hom_st)

print(f"Shape of {cond1} anndata: ")
print(data_wt_st.shape)


print(f"Shape of {cond2} anndata: ")
print(data_hom_st.shape)


#do wilcoxon Rank Sum test for each gene in 14832
print("Stats testing for each gene")
genes_stat = []
p_vals_stat = []
for i in range(len(data_wt_st.var_names)):
    if i % 1000 == 0:
    	print(".", end="")
    mwu = sp.stats.mannwhitneyu(data_wt_st.X[:,i], data_hom_st.X[:,i], alternative='two-sided')
    p_vals_stat.append(mwu[1])
    genes_stat.append(data_wt_st.var_names[i])
print(".")
adj_p_vals_stat = mt.multipletests(pvals=p_vals_stat, alpha=0.05, method='fdr_bh')[1]

#ma plots
m = []
a = []

m_sig = []
a_sig = []

for i in range(len(data_wt_st.var_names)):
    mean_wt = np.mean(data_wt_st.X[:,i])
    mean_hom = np.mean(data_hom_st.X[:,i])
    mean = np.log2(mean_hom/mean_wt)
    #av = mean_wt*mean_hom 
    av = mean_wt
    m.append(mean)
    a.append(av)
    
    if adj_p_vals_stat[i] < 0.00000000000000000005:
        m_sig.append(mean)
        a_sig.append(av)

# Create dataframe
df_stat = pd.DataFrame(index=data_wt_st.var_names)
df_stat['name'] = data_wt_st.var_names
df_stat['p_value'] = p_vals_stat
df_stat['adj_p_value'] = adj_p_vals_stat
m = np.array(m)
df_stat['log2_FC'] = m
updown = ['Up' if m[i] > 0 else 'Down' for i in range(len(m))]
df_stat['direction'] = updown
df_stat = df_stat.sort_values('adj_p_value')

# Pickling them
bn = cond1 + "_" + cond2 + "_" + clus 
df_stat.to_csv(bn + '.dvg.csv')
df_stat.to_pickle(bn + '.dvg.pkl')