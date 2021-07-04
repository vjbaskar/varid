#!/usr/bin/env python3


"""data
Usage: 04.plot_expAndVar.py [--ep pickefile1] [--vp pickefile2] [--c1 condition1] [--c2 condition2] [-s subset_name] [-g genefile ] [-h help]


Options:
    -h help
    --c1=condition 1 name ideally wildtype
    --c2=condition 2 name ideally mutant
    -s=cluster/group of cells to subset (eg. HSCs)
    -g=genefile containing one gene per row (eg. genelist1.txt)
    --ep=pickefile for expression deg
    --vp=picklefile for variability dvg
"""



from docopt import docopt


genefile = 'ery_GO_0030218.txt'
cond1 = 'WT'
cond2 = 'Jak2Hom'
clus = 'HSCs'
deg_pkl = '../varid_1/' + cond1 + '_' + cond2 + '_' + clus + '.deg.pkl'
dvg_pkl = '../varid_1/' + cond1 + '_' + cond2 + '_' + clus + '.dvg.pkl'

arguments = docopt(__doc__)

cond1 = arguments['--c1']
cond2 = arguments['--c2']
clus = arguments['-s']
genefile = arguments['-g']
deg_pkl = arguments['--ep']
dvg_pkl = arguments['--vp']

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
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list(name='gene_cmap',
                                        colors=[(0, 'lightgrey'), (0.1, 'yellow'), (0.5, 'red'), (1, 'darkred')])


def filter_df(df, genes, exp_lfc = 0.263, exp_pval = 0.05, var_lfc = 0.263, var_pval = 0.05):
    print(len(genes))
    df_genes = df[df.index.isin(genes)]
    df_filt = df_genes[df_genes['adj_p_value'] <= var_pval ]
    #df_filt = df_filt[np.abs(df_filt['log2_FC']) >= var_lfc ] 
    #df_filt = df_filt[np.abs(df_filt['logfoldchanges']) >= exp_lfc ]
    #df_filt = df_filt[df_filt['pvals_adj'] <= exp_pval ] 
    #df_filt = df_filt[df_filt['pvals_adj'] <= 0.05 ]
    #print(df_filt.shape)
    return df_filt
def ci(arr, inter):
    tmp = st.t.interval(inter, len(arr)-1, loc=np.mean(arr), scale=st.sem(arr))
    return (tmp[1]-tmp[0])/2

def plotting_devg(df_filt):
    var_mean = np.mean(df_filt['log2_FC'])
    expr_mean = np.mean(df_filt['logfoldchanges'])
    var_ci = ci(df_filt['log2_FC'],.95 )
    expr_ci = ci(df_filt['logfoldchanges'],.95 )
    
    plt.figure(figsize=(8,5))
    x = 0
    plt.errorbar(x+0.1, var_mean, yerr=var_ci, fmt='.', color='purple', capsize=5, label='Variability', markersize=20, elinewidth=3, capthick=3)
    plt.errorbar(x-0.1, expr_mean, yerr=expr_ci, fmt='.', color='blue', capsize=5, label='Expression', markersize=20, elinewidth=3, capthick=3)
    plt.hlines(0,-.5,1.5, linestyle='dashed')
    plt.xticks(ticks=[0],labels=[cond1 + "_" + cond2], rotation=90, fontsize=16)
    #plt.xlim(-.5,7.5)
    plt.legend(fontsize=16)
    #plt.savefig('/home/sw631/thesis_figs/chap4/varID/pointplot_ery_genes.svg', dpi=300)
    #plt.show()
    return plt, var_mean, expr_mean, var_ci, expr_ci

with open(genefile, 'r') as f:
    genes = [ i.strip() for i in f.readlines() ]
print('Total genes in file = {}'.format(len(genes)))

with open(deg_pkl, 'rb') as f:
    de = pickle.load(f)
print("{} shape is {}".format(cond1, de.shape))
de.index = list(de['names'])

with open(dvg_pkl, 'rb') as f:
    df = pickle.load(f)
print("{} shape is {}".format(cond2, df.shape))

de_df = de.merge(df, left_index=True, right_index=True)
de_df = de_df.sort_values('adj_p_value')
print("merged df shape is {}".format(de_df.shape))

df_filt = filter_df(de_df, genes) #, exp_lfc = 0.1, var_lfc = 0.263)
print("After filtering df shape is {}".format(df_filt.shape))
total_genes = df_filt.shape[0]


ofile = cond1+"_"+cond2+"_"+ "clus" 
with PdfPages(ofile + "_evar.pdf") as pdf:
    plt, var_mean, expr_mean, var_ci, expr_ci = plotting_devg(df_filt)
    plt.title(cond1 + " " + cond2 + "| Total genes = " + str(total_genes))
    pdf.savefig()

df = pd.DataFrame({"c1": cond1, "c2": cond2, "clus": clus, 
                   "var_mean": var_mean, "expr_mean": expr_mean, 
                   "var_ci": var_ci, "expr_ci": expr_ci, 'geneset': genefile, "total_genes": len(genes), 'genes_plotted': total_genes}, index = [ofile ])

df.to_csv(ofile + '_evars.csv')




