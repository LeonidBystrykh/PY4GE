#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 10:27:08 2020
try second example, Kowalczyk 3 cell types
@author: leonid
Useful options in
https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scanpy/scanpy_05_dge.html
"""

import scanpy as sc, pandas as pd, seaborn as sns, matplotlib.pyplot as plt, numpy as np
#%%
#df=pd.read_csv("GSE61533_HTSEQ_counts.txt", sep="\t")
# scanpy has own function to open the file
dir="/home/leonid/Documents/review_SC_python/"
adata=sc.read_csv(dir + "GSE59114_C57BL6_Y_all.csv", delimiter="\t").T

#%%
adata.obs
#%%
adata.var
#%%
adata.X

#%%
plt.hist(adata.X[0])
#%%
data=adata[:,'Cd34']
print(data.X)
#%%
CD_genes = list(adata.var_names.str.startswith('Cd'))
print(adata[:,CD_genes].var)
#%%
#plot the hiest expression, note the table orientation
sc.pl.highest_expr_genes(adata, n_top=20, )
#%% restore read counts
adata.X=2**adata.X-1
#check mitochondrial genes, HTSEQ
mito_genes = adata.var_names.str.startswith('mt-')
mito_genes
print(list(mito_genes).count(True))
#estimate counts per cell
adata.obs['n_counts'] = adata.X.sum(axis=1)
#looks like filtering
sc.pp.filter_cells(adata, min_genes=1000)
sc.pp.filter_genes(adata, min_cells=5)
#another filtering
#adata = adata[adata.obs.n_genes < 2500, :]
#visuals
#%%
sc.pl.scatter(adata, x='n_counts', y='n_genes')
sc.pl.violin(adata, ['n_genes','n_counts'], jitter=0.4, multi_panel=True)
#%%
#normalise and log-transform
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

#%%
#find and plot variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.3, min_disp=0.3)
sc.pl.highly_variable_genes(adata)
#print highly variable genes
adata = adata[:, adata.var.highly_variable]
adata.var.highly_variable[:50]
#%%
print('Variable genes total ', len(adata.var.highly_variable))
#%%
sc.pp.scale(adata, max_value=1) #in seurat it comes after findVariableFeatures
#%%
#need for DE
adata.raw = adata
#%%
#PCA and attempt to cluster, color by gene expressions
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color="Npl")
#%%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['Flt3'])

#%% run leiden clustering
sc.tl.leiden(adata)
#%%
sc.pl.umap(adata, color='leiden')
#%% make groups by cell type
cell_types=[]
for i in range(len(adata.obs['leiden'])):
    print(adata.obs_names[i], adata.obs['leiden'][i])
    if 'LT' in adata.obs_names[i]:
        cell_types.append(0)
    if 'ST' in adata.obs_names[i]:
        cell_types.append(1)  
    if 'MPP' in adata.obs_names[i]:
        cell_types.append(2)
#%% add cell_type categories
adata.obs['cell_types']=cell_types
adata.obs['cell_types']=adata.obs['cell_types'].astype('category')
#%%
#sc.pl.umap(adata, color=['cell_types','leiden'])
sc.pl.umap(adata, color=['cell_types'])
#%%
sc.pl.pca(adata, color=['cell_types','leiden']) 
#%%
sc.pl.pca(adata, color=['cell_types'], components='1,2')
#%%
#out=pd.DataFrame([adata.uns['pca'][0],adata.uns['pca'][1],cell_types])
#%% you can use either leiden clusters or cell_types
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
#%%
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
#%% more specific way of finding markers
#LT: 1,5
#ST: 0,4
#MPP: 2,3
#sc.tl.rank_genes_groups(adata, 'leiden', groups=['1','4'], 
 #                       reference='0', 
  #                      method='wilcoxon')
#sc.pl.rank_genes_groups(adata, groups=['1','4'], n_genes=20)
#%% one way of finding all markers
sc.tl.rank_genes_groups(adata, 'cell_types', method='wilcoxon')
#%%
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
#%%
result=adata.uns['rank_genes_groups']['names']
pvals=adata.uns['rank_genes_groups']['pvals_adj']
for i in range(40):
    print(result[i], pvals[i])
#it is more or less all you can do with results
#there are more fancy options for the illustrations, check ref above
#%% other way of finding all markers, but two methods do not match
    # this option hase lower match to Seurat list
LT_DE = sc.get.rank_genes_groups_df(adata, group=['0','1','2'],# method='wilcoxon', 
                                    pval_cutoff=0.01, log2fc_min=1)#['names']
#%%
LT_DE.to_csv(dir+ "Kow_scanpy_DE_0303.tsv", sep="\t")