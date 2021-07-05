# Differential variability calculator

## Requirements
* scanpy
* numpy
* scipy
* anndata
* matplotlib
* pickle
* scikit-learn
* docopt
* pandas

## Example
Let us assume that you have a scanpy anndata `cluster_nikiLabelled.h5ad` and there are two conditions `WT` and `Jak2Hom`. 

### Step 1: Quantile norm variabilities
Generate quantile normalised variabilities for all the genes across all cells. 

```
python 01.compute_varid.py -d cluster_nikiLabelled.h5ad -a WT -b Jak2Hom
mkdir WT_Jak2Hom
mv *corrections.pdf *var.h5ad *qnorm.h5ad WT_Jak2Hom
```

### Step 2: Diff Variable Genes (dvgs)
For conditions WT and Jak2Hom there is a cluster called `HSCs` in `WT_qnorm.obs['niki']` and in `Jak2Hom_qnorm.obs['niki']`.
If you want to consider all cells, then artificially add an obs column with the same value and use that.

```
python 02.dvg.py --a1 WT_Jak2Hom/WT_qnorm.h5ad --a2 WT_Jak2Hom/Jak2Hom_qnorm.h5ad --c1 WT --c2 Jak2Hom -o niki -s HSCs
```

### Step 3: Diff Expressed Genes (degs)
Ideally Step 2 marks the completion of everything related to variability analysis. However, in order to compare degs with dvgs this script is useful. This calculates differentially expressed genes.

```
python 03.deg.py -d cluster_nikiLabelled.h5ad --c1 WT --c2 Tet2Hom -s HSCs -o niki
```

### Step 4: 

Just a small plotting script for plotting degs and dvgs for a given gene set `ery_GO_0030218.txt`.

```
python 04.plot_expAndVar.py --c1 WT --c2 Jak2Hom -s HSCs -g ery_GO_0030218.txt --ep WT_Jak2Hom_HSCs.deg.pkl --vp WT_Jak2Hom_HSCs.dvg.pkl
```

#### Docker container
If you want to use docker container try `docker pull vjbaskar/varid`. Or build in singularity using the command below.
```
singularity build varid.sif docker://vjbaskar/varid
```
All the python scripts are present in the `/` folder. For example, 
```
singularity run varid.sif python /04.plot_expAndVar.py --c1 WT --c2 Jak2Hom -s HSCs -g ery_GO_0030218.txt --ep WT_Jak2Hom_HSCs.deg.pkl --vp WT_Jak2Hom_HSCs.dvg.pkl
```
You can use all the above steps with just substituting `*.py` to `/*.py`
