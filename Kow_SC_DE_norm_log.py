#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 07:30:07 2020
try filter Wilson data pandas style
do first normalise then log (as it is done in seurat and scanpy)
@author: leonid
"""
print("import packages")
import pandas as pd, numpy as np, matplotlib.pyplot as plt, scipy.stats as ss, math
import scipy.spatial.distance as distance
from sklearn.preprocessing import normalize
from umap import UMAP
from sklearn.manifold import TSNE, SpectralEmbedding
#import seaborn as sns
#%%
my_d="/home/leonid/Documents/review_SC_python/"
df_SC=pd.read_csv(my_d + "GSE59114_C57BL6_Y_all.csv", sep="\t", index_col="Gene Symbol")
#convert log-transform data back to reads
new=(2**df_SC.values-1)
# restore data frame
df_SC=pd.DataFrame(new, index=df_SC.index, columns=df_SC.columns)
#%%
df_SC.head()#get some information about our Data-Set
#%%
df_SC.describe()#stat summary
#%% change column names, collect cell type data
new,groups=[],[]
for item in df_SC.columns:
    new.append(item.replace('young_',''))
    if 'LT' in item:
        groups.append(0)
    if 'ST' in item:
        groups.append(1)
    if 'MPP' in item:
        groups.append(2)
df_SC.columns=new
#%% count non-zero values
nonz_col=df_SC.astype(bool).sum(axis=0)
nonz_row=df_SC.astype(bool).sum(axis=1)
to_remove=[x for x in nonz_row if x<5]
print("rows to remove", len(to_remove))
to_remove=[x for x in nonz_col if x<2000]
print("cols to remove", len(to_remove))
#%%
plt.subplot(121)
plt.hist(nonz_col)
plt.title("per columns")
plt.subplot(122)
plt.hist(nonz_row)
plt.title("per rows")
#%% remove genes with low presense (in rows)
df_SC["present"]=nonz_row
df_SC=df_SC[df_SC['present']>5]
df_SC=df_SC.drop(labels='present',axis=1)
#%%
totals=df_SC.sum(axis=0)
present=df_SC.astype(bool).sum(axis=0) #this is present in columns
#%%
plt.scatter(totals,present, c=groups)
plt.xlabel('totals')
plt.ylabel('present')
#%%
# we could try some normalizations from sklearn, but RNAseq specific case is not there
#besides, we better preserve df format
filtered=df_SC #rename in case you need SC again
medianT=np.median(totals)
medianP=np.median(present)
for i in range(0,len(totals)):
    data=filtered[filtered.columns[i]]
    step1=data*medianT/totals[i]
    step2=step1*present[i]/medianP
    step3=np.log2(step2+1)
    filtered[filtered.columns[i]]=step3
#%% alternative normalization with sklearn normalize
    #with the same filter thresholds makes more genes on the list
    #it makes results less similar to Seurat, missing Apoe, Mpo
    #probably less preferrable
#norm2=normalize(df_SC, norm='l2', axis=0)*10000+1 #0 is for columns
#norm_df=pd.DataFrame(norm2, columns=df_SC.columns, index=df_SC.index)
#filtered=norm_df.applymap(np.log2)
#%% Now we are ready for the next step
#try to visualise variations per gene and pinpoint the most variable
filtered["Vars"]=filtered[filtered.columns].var(axis=1)
plt.hist(filtered["Vars"])
#%% simple way of filtering by value
filtered=filtered.loc[filtered['Vars'] > 3] #threshold is based on the histogram data

#%% create means column
filtered['means'] = filtered[filtered.columns].mean(axis=1)
plt.hist(filtered['means'])
#%% filter by means
filtered=filtered.loc[filtered['means'] > 3] #threshold is based on the histogram data
#%% remove Vars and means columns
new=filtered.drop(labels=['Vars','means'],axis=1)
#%% scale gene expressions
scaled=normalize(new, norm='l2', axis=1) #1 is for rows
df=pd.DataFrame(scaled, index=new.index, columns=new.columns)
#%%
# alternatively we can use sklearn filter low variance
#however it gives the same result and take extra effort to preserve all names
#it naturally returns np.array type
#import sklearn.feature_selection as skf
# function from https://stackoverflow.com/questions/39812885/retain-feature-names-after-scikit-feature-selection
#def variance_threshold_selector(data, threshold):
#    selector = skf.VarianceThreshold(threshold)
#    selector.fit(data)
#    return data[data.columns[selector.get_support(indices=True)]]
# run the filter by variance
#new=variance_threshold_selector(logged.T, 3).T

#%%
#try PCA
cells=df.columns
X=scaled.T
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pComp = pca.fit_transform(X).T
print(pComp)
print("Explained variance")
print(pca.explained_variance_ratio_)
plt.scatter(pComp[0],pComp[1], c=groups)
#%% save the PCA data
my_array=np.array([groups,pComp[0],pComp[1]]).T
out=pd.DataFrame(my_array, columns=["groups", "PC1", "PC2"])
out.to_csv(my_d+ "Kow_pd_sk_norm_PCA.tsv", sep="\t")
#%% add text if you want
for i in range(len(cells))[:10]:
    plt.text(pComp[0][i],pComp[1][i], cells[i])
#%%
#try tSNE

'''
here we can play with several parameters, such as
perplexity: default 30, less is more variable, more is less
metric: Valid values for metric are:
From scikit-learn: [‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’]. 
These metrics support sparse matrix inputs.
From scipy.spatial.distance: [‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’, ‘hamming’, 
‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, 
‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’] 
See the documentation for scipy.spatial.distance for details on these metrics.
'''
from sklearn.metrics import pairwise_distances as pwd
dist=pwd(X, Y=None, metric='correlation')
#%%
tsne = TSNE(n_components=2, #init='pca', 
            learning_rate=100,
            verbose=5,
            perplexity=15, metric='precomputed', method='exact')
trans_data = tsne.fit_transform(dist).T
print("t-SNE: %.2g sec")
plt.scatter(trans_data[0], trans_data[1], c=groups)
plt.title("t-SNE")
#%%
# affinity{‘nearest_neighbors’, ‘rbf’, ‘precomputed’, 
embedding = SpectralEmbedding(n_components=2, affinity='rbf')
X_transformed = embedding.fit_transform(X).T
plt.scatter(X_transformed[0], X_transformed[1], c=groups)
plt.title("Spectral Embedding")
#%%
reducer = UMAP(n_neighbors=15, min_dist=0.08,local_connectivity=10)
embedding = reducer.fit_transform(X)
#embedding.shape
plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=groups,
   # c=clf.predict(X)#colorz clf is classifier in random forest, use it after
    ) 
plt.gca().set_aspect('equal', 'datalim')
plt.xlabel("UMAP_1")
plt.ylabel("UMAP_2")
#plt.title('UMAP Kowalczyk data')
#%%
from sklearn.cluster import DBSCAN
clustering = DBSCAN(eps=3, min_samples=2).fit(embedding)
colorz=list(clustering.labels_)
for i in range(len(colorz)):
    plt.text(embedding[:, 0][i],
    embedding[:, 1][i], colorz[i])
#%% define groups and try stats for DE
LT=[x for x in df.columns if "LT" in x]
ST=[x for x in df.columns if "ST" in x]
MPP=[x for x in df.columns if "MPP" in x]
#%%
counter=0
for i in range(len(df.index)):
    group1=df[LT]
    group2=df[ST]
    group3=df[MPP]
  #  stat, pval=ss.ttest_ind(group1.loc[df.index[i]],group2.loc[df.index[i]]) #for two groups only
   # stat, pval=ss.f_oneway(group1.loc[df.index[i]],group2.loc[df.index[i]],group3.loc[df.index[i]])
    stat, pval=ss.kruskal(group1.loc[df.index[i]],group2.loc[df.index[i]],group3.loc[df.index[i]])
    if pval<1e-12:
        counter+=1
        print(counter,df.index[i], pval)
#%%
#try random forest to fit DBSCAN and find genes
from sklearn.ensemble import RandomForestClassifier
y=groups
clf = RandomForestClassifier(max_depth=2, random_state=0)
clf.fit(X, y)
#%%
print(clf.predict(X)) #predict all data
print("success score")
print(clf.score(X,y))
imps = clf.feature_importances_
print("important features")
for i in range(len(imps)):
    if imps[i]>0.005:
        print(i,imps[i], df.index[i])
#%%
#and we can select other groups by high-low expression of some genes or groups of genes
        #find differences between those groups and so on.
        #it will be very much the same kind of lines as already shown above
        