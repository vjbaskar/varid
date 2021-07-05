#!/usr/bin/env python3


"""data
Usage: 01.compute_varid.py [-d datafile] [-a condition1] [-b condition2] [-h help]


Options:
    -h --help
    -d=<datafiile>  anndatafile
    -a=<condition1> condition 1
    -b=<condition2> condition 2

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
import scipy
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list(name='gene_cmap',
                                        colors=[(0, 'lightgrey'), (0.1, 'yellow'), (0.5, 'red'), (1, 'darkred')])

anndatafile = '../ANALYSIS1/cluster_nikiLabelled.h5ad'
cond1 = 'WT'
cond2 = 'Tet2Hom'

arguments = docopt(__doc__)
anndatafile = arguments["-d"]
cond1 = arguments['-a']
cond2 = arguments['-b']


def filter_genes(adata, x):
    boo = np.sum(adata.X > x, axis=0) > x
    print('Total number of genes passing filtering is:', sum(boo))
    adata = adata[:, boo]
    return adata

def calculate_correction(adata, d=2):
    m = np.mean(adata.X, axis=0)
    v = np.var(adata.X, axis=0)
    m = np.log2(m)
    v = np.log2(v)
    r = np.polyfit(m, v, d)
    return r

def plot_correction(adata, r):
    plt.figure(figsize=(5, 5))
    m = np.mean(adata.X, axis=0)
    v = np.var(adata.X, axis=0)
    m = np.log2(m)
    v = np.log2(v)
    x = np.linspace(-8, 8, 100)
    y = r[0]*x**2 + r[1]*x + r[2]
    plt.scatter(m, v, alpha=0.1)
    plt.scatter(x, y)
    #plt.show()
    return plt

def calculate_variabilities(adata, D, r, k=10):
    R = np.argpartition(D, k, axis=1)[:,:k]
    variability = np.zeros(adata.shape)
    for i in range(adata.shape[0]):
        lnexpr = adata.X[R[i],:]
        lvars_unc = np.var(lnexpr, axis=0)
        lmeans= np.mean(lnexpr, axis=0) + 0.05
        varfac = 2**(r[0]*np.log2(lmeans)**2 + r[1]*np.log2(lmeans) + r[2])
        lvars = lvars_unc/varfac
        variability[i,:] = lvars
    return variability


def retain_shared_genes(data1, data2):
    boo1 = [name in data2.var_names for name in data1.var_names]
    boo2 = [name in data1.var_names for name in data2.var_names]
    data1 = data1[:,boo1]
    data2 = data2[:,boo2]
    print(len(data1.var_names))
    print(len(data2.var_names))
    print(data1.var_names == data2.var_names)
    return data1, data2

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

# Read anndata
print("Reading anndata")
adata = sc.read(anndatafile)
sc.pp.normalize_total(adata, exclude_highly_expressed=True, target_sum=10000)

if scipy.sparse.issparse(adata.X):
    print("Data being converted to dense")
    adata.X = adata.X.todense()
data_wt = adata[adata.obs['condition']==cond1,:].copy()
data_hom = adata[adata.obs['condition']==cond2,:].copy()

#data_wt.X = data_wt.X.A
#data_hom.X = data_hom.X.A


# Filter genes
print("Filtering low information genes")
data_wt = filter_genes(data_wt, 3)
data_hom = filter_genes(data_hom, 3)


# Mean-variance fitting
print("Fitting mean-variance")
r_wt = calculate_correction(data_wt)
r_hom = calculate_correction(data_hom)

with PdfPages(cond1+"_"+cond2+"_corrections.pdf") as pdf:
    plt_wt = plot_correction(data_wt, r_wt)
    plt_wt.title(cond1)
    pdf.savefig()
    plt_wt.close()
    plt_hom = plot_correction(data_hom, r_hom)
    plt_hom.title(cond2)
    pdf.savefig()
    plt_hom.close()


data_wt = anndata.AnnData(X=data_wt.X, obs=data_wt.obs, var=data_wt.var)
data_hom = anndata.AnnData(X=data_hom.X, obs=data_hom.obs, var=data_hom.var)

print(f"Calculating PW distances for {cond1}")
sc.pp.pca(data_wt, n_comps=50)
X_pca = data_wt.obsm['X_pca']
D_wt = sp.spatial.distance.squareform(sp.spatial.distance.pdist(X_pca, metric='euclidean'))


print(f"Calculating PW distances for {cond2}")
sc.pp.pca(data_hom, n_comps=50)
X_pca = data_hom.obsm['X_pca']
D_hom = sp.spatial.distance.squareform(sp.spatial.distance.pdist(X_pca, metric='euclidean'))


print(f"Calculating variability matrix for {cond1}")
var_wt = calculate_variabilities(data_wt, D_wt, r_wt)
print(f"Calculating variability matrix for {cond2}")
var_hom = calculate_variabilities(data_hom, D_hom, r_hom)

data_wt_var = anndata.AnnData(X=var_wt, obs=data_wt.obs, var=data_wt.var)
data_hom_var = anndata.AnnData(X=var_hom, obs=data_hom.obs, var=data_hom.var)

print("Saving anndata")
data_wt_var.write(cond1 + '_var.h5ad')
data_hom_var.write(cond2 + '_var.h5ad')

#### Quantile normalisation

data_wt = data_wt_var.copy()
data_hom = data_hom_var.copy()

print("Retaining only shared genes")
#subset genes in each experiment to only shared genes
data_wt, data_hom = retain_shared_genes(data_wt, data_hom)

print("Performing quantile normalisation")
#Quantile Normalise?
data = np.vstack((data_wt.X,data_hom.X)).T
ranks = np.zeros(data.shape)
sorts = np.zeros(data.shape)
qn_dat = np.zeros(data.shape)
for i in range(data.shape[1]):
    #print(i)
    rankswt = data[:,i].argsort().argsort()
    sortwt = np.array(sorted(data[:,i]))
    ranks[:,i] = rankswt
    sorts[:,i] = sortwt

means = np.mean(sorts, axis=1)
for i in range(data.shape[1]):
    #print(i)
    qn_dat[:,i] = means[np.array([int(ranks[:,i][j]) for j in range(data.shape[0])])]
    
data_wt = anndata.AnnData(X=qn_dat[:,:data_wt.shape[0]].T, obs=data_wt.obs, var=data_wt.var)
data_hom = anndata.AnnData(X=qn_dat[:,data_wt.shape[0]:].T, obs=data_hom.obs, var=data_hom.var)

print("Saving anndata")
data_wt.write(cond1 + '_qnorm.h5ad')
data_hom.write(cond2 + '_qnorm.h5ad')
