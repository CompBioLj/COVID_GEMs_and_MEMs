# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 13:31:11 2021

@author: alex
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr, spearmanr
from itertools import permutations, product, combinations


d=pd.read_csv('M_reaction_end_all.csv',sep=',')

col=d.columns.tolist()
col.pop(0)
M=d[col]

M=np.transpose(M.values)

uninfected =col[:20]
infected=col[20:]
infection = (uninfected, infected)

method1 =['Tinit_NHBE_H', 'Tinit_NHBE_I','Tinit_A549_H','Tinit_A549_I','Tinit_293T_H','Tinit_293T_I','Tinit_Lung_H','Tinit_Lung_I','Tinit_CALU_H','Tinit_CALU_I']
method2 =['Gimme_NHBE_H', 'Gimme_NHBE_I','Gimme_A549_H','Gimme_A549_I','Gimme_293T_H','Gimme_293T_I','Gimme_Lung_H','Gimme_Lung_I','Gimme_CALU_H','Gimme_CALU_I']
method3 =['iMAT_NHBE_H', 'iMAT_NHBE_I','iMAT_A549_H','iMAT_A549_I','iMAT_293T_H','iMAT_293T_I','iMAT_Lung_H','iMAT_Lung_I','iMAT_CALU_H','iMAT_CALU_I']
method4 =['init_NHBE_H', 'init_NHBE_I','init_A549_H','init_A549_I','init_293T_H','init_293T_I','init_Lung_H','init_Lung_I','init_CALU_H','init_CALU_I']

MEM = (method1, method2, method3, method4)

Dat1=['Tinit_NHBE_H', 'Tinit_NHBE_I','init_NHBE_H', 'init_NHBE_I','Gimme_NHBE_H', 'Gimme_NHBE_I','iMAT_NHBE_H', 'iMAT_NHBE_I']
Dat2=['Tinit_CALU_H', 'Tinit_CALU_I','init_CALU_H', 'init_CALU_I','Gimme_CALU_H', 'Gimme_CALU_I','iMAT_CALU_H', 'iMAT_CALU_I']
Dat3=['Tinit_A549_H', 'Tinit_A549_I','init_A549_H', 'init_A549_I','Gimme_A549_H', 'Gimme_A549_I','iMAT_A549_H', 'iMAT_A549_I']
Dat4=['Tinit_293T_H', 'Tinit_293T_I','init_293T_H', 'init_293T_I','Gimme_293T_H', 'Gimme_293T_I','iMAT_293T_H', 'iMAT_293T_I']
Dat5=['Tinit_Lung_H', 'Tinit_Lung_I','init_Lung_H', 'init_Lung_I','Gimme_Lung_H', 'Gimme_Lung_I','iMAT_Lung_H', 'iMAT_Lung_I']

Data=(Dat1,Dat2,Dat3,Dat4,Dat5)

groups = {"infection": infection, "MEM": MEM, "Data":Data}
labels = {"infection": ("uninfected","infected"),
          "MEM" : ("Tinit", "Gimme", "iMAT", "init"),
          "Data": ("NHBE","CALU","A549","293T","Lung")}



M = M[:,~np.all(M==0, axis=0)]
M = M[:,~np.all(M==1, axis=0)]
M = np.array(M, dtype=float)

M = StandardScaler(with_mean=True, with_std=False).fit_transform(M) #centering only

models=list(col)

#pca = PCA(n_components=3)
pca = PCA(n_components=2)
comps = pca.fit_transform(M)
#df= pd.DataFrame(data = comps, columns = ['PC1', 'PC2', 'PC3'])
df= pd.DataFrame(data = comps, columns = ['PC1', 'PC2'])
df['models'] = models
df_pca = pd.DataFrame()
df_pca = df_pca.append(pd.DataFrame(pca.explained_variance_).T)
df_pca = df_pca.append(pd.DataFrame(pca.explained_variance_ratio_).T)
#df_pca.columns = ['PC1', 'PC2', 'PC3']
df_pca.columns = ['PC1', 'PC2']
df_pca['label'] = ['explained variance', 'explained variance ratio']
df_pca = df_pca.set_index('label')
#df_pca.to_csv("explained_variance.csv")
pca_explained = pca.explained_variance_ratio_
#for c in combinations(range(3), 2):
for c in combinations(range(2), 2):
    i1 = c[0]
    i2 = c[1]
    for group in groups:        
        for subgroup, label in zip(groups[group], labels[group]):
            locs = np.isin(models, subgroup)
            plt.plot(comps[locs,i1], comps[locs,i2],"o", label=label)

        """
        for model, x, y in zip(models, comps[:,i1], comps[:,i2]):
            #plt.text(x,y,model[6:])
            for subgroup in groups[group]:
                if model in subgroup:            
                    plt.text(x+0.05,y+0.05,model)
        """

        plt.title(group )
        plt.xlabel("PC"+str(i1+1) + " (" + str(round(100*pca_explained[i1],2))+"%)")
        plt.ylabel("PC"+str(i2+1) + " (" + str(round(100*pca_explained[i2],2))+"%)")
        plt.legend()
        plt.gcf().set_size_inches(20,10)
        #plt.savefig("_PC"+str(i1+1)+'_'+"PC"+str(i2+1)+'_'+group+".pdf", format="pdf", bbox_inches = 'tight')
        plt.savefig("figure/All/PCA/ALL_PC"+str(i1+1)+'_'+"PC"+str(i2+1)+'_'+group+".pdf", format="pdf", bbox_inches = 'tight')
        plt.show()
        
factors = list(groups.keys())
#Rs = np.zeros((len(factors), 3))
#rhos = np.zeros((len(factors), 3))
Rs = np.zeros((len(factors), 2))
rhos = np.zeros((len(factors), 2))


for ii, factor in enumerate(groups):
    scores1 = []
    scores2 = []
    #scores3 = []
    for i in range(len(groups[factor])):
        idxs = np.array(np.where(np.isin(models, groups[factor][i])==True)).flatten()
        scores1.append(sorted(df.iloc[idxs, 0].values))
        scores2.append(sorted(df.iloc[idxs, 1].values))
        #scores3.append(sorted(df.iloc[idxs, 2].values))
        
    for idx in permutations(range(len(scores1))):
        s1 = []
        s2 = []
        #s3 = []        
        for i in idx:
            s1 += scores1[i]
            s2 += scores2[i]
            #s3 += scores3[i]
            
        
        R_PC1 = pearsonr(np.arange(len(s1)), s1)[0]
        R_PC2 = pearsonr(np.arange(len(s2)), s2)[0]
        #R_PC3 = pearsonr(np.arange(len(s3)), s3)[0]

        rho_PC1 = spearmanr(np.arange(len(s1)), s1)[0]
        rho_PC2 = spearmanr(np.arange(len(s2)), s2)[0]
        #rho_PC3 = spearmanr(np.arange(len(s3)), s3)[0]

        Rs[ii, 0] = max(Rs[ii, 0], abs(R_PC1))
        Rs[ii, 1] = max(Rs[ii, 1], abs(R_PC2))
        #Rs[ii, 2] = max(Rs[ii, 2], abs(R_PC3))


Rs = Rs ** 2

#df_R = pd.DataFrame(columns = ["PC1", "PC2", "PC3"], data=Rs)
df_R = pd.DataFrame(columns = ["PC1", "PC2"], data=Rs)

df_R['factor'] = factors
df_R = df_R.set_index("factor")
df_R.to_csv("figure/CSV/explained_variance_R.csv")

# tsne
n_components= 2
perplexity = 15 # change the perplexity 5 to 20
tSNE = TSNE(n_components=n_components, perplexity = perplexity,  n_iter=5000)
comps = tSNE.fit_transform(M)
df = pd.DataFrame(data = comps, columns = ['PC1', 'PC2', 'PC3'][:n_components])
df['model'] = models

for c in combinations(range(n_components), 2):
    i1 = c[0]
    i2 = c[1]
    for group in groups:
        for subgroup, label in zip(groups[group], labels[group]):
            locs = np.isin(models, subgroup)
            plt.plot(comps[locs,i1], comps[locs,i2],"o", label=label)
            
        #for model, x, y in zip(models, comps[:,i1], comps[:,i2]):        
        #    plt.text(x+0.05,y+0.05,model)

        plt.title(group)
        plt.xlabel("C"+str(i1+1))#plt.xlabel("PC"+str(i1+1) + " (" + str(round(100*pca_explained[i1],2))+"%)")
        plt.ylabel("C"+str(i2+1))#plt.ylabel("PC"+str(i2+1) + " (" + str(round(100*pca_explained[i2],2))+"%)")
        plt.legend()
        plt.gcf().set_size_inches(20,10)
        #plt.savefig("C"+str(i1+1)+'_'+"C"+str(i2+1)+'_'+group+".pdf", format="pdf", bbox_inches = 'tight')
        #plt.savefig("figure/All/TSNE/All_Tsne_C"+str(i1+1)+'_'+"C"+str(i2+1)+'_'+group+".png", format="png", bbox_inches = 'tight')
        plt.savefig("figure/All/TSNE/All_Tsne_C"+str(i1+1)+'_'+"C"+str(i2+1)+'_'+group+".pdf", format="pdf", bbox_inches = 'tight')
        plt.show()
