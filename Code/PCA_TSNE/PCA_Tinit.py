# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 13:33:14 2021

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

#Tinit
d=pd.read_csv('M_reaction_end_Tinit.csv',sep=',')

col_Tinit=['Tinit_NHBE_H','Tinit_A549_H', 'Tinit_293T_H','Tinit_Lung_H','Tinit_CALU_H','Tinit_NHBE_I','Tinit_A549_I','Tinit_293T_I','Tinit_Lung_I','Tinit_CALU_I']
M_Tinit=d[col_Tinit]
M_Tinit=np.transpose(M_Tinit.values)



data_Tinit_NHBE=['Tinit_NHBE_H','Tinit_NHBE_I']
data_Tinit_A549=['Tinit_A549_H','Tinit_A549_I']
data_Tinit_293T=['Tinit_293T_H','Tinit_293T_I']
data_Tinit_CALU=['Tinit_CALU_H','Tinit_CALU_I']
data_Tinit_Lung=['Tinit_Lung_H','Tinit_Lung_I']
data_Tinit=(data_Tinit_NHBE,data_Tinit_A549,data_Tinit_293T,data_Tinit_CALU,data_Tinit_Lung)


uninfected =col_Tinit[:5]
infected=col_Tinit[5:]
infection = (uninfected, infected)

groups = {"infection": infection, "Data":data_Tinit}
labels = {"infection": ("uninfected","infected"),
           "Data": ("NHBE","CALU","A549","293T","Lung")}


M_Tinit = M_Tinit[:,~np.all(M_Tinit==0, axis=0)]
M_Tinit= M_Tinit[:,~np.all(M_Tinit==1, axis=0)]
M_Tinit = np.array(M_Tinit, dtype=float)

M_Tinit = StandardScaler(with_mean=True, with_std=False).fit_transform(M_Tinit) #centering only

models=list(col_Tinit)

#pca = PCA(n_components=3)
pca = PCA(n_components=2)
comps = pca.fit_transform(M_Tinit)
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
        plt.savefig("figure/tINIT_plot/tINIT_PC"+str(i1+1)+'_'+"PC"+str(i2+1)+'_'+group+".pdf",bbox_inches = 'tight')
        plt.show()
        
factors = list(groups.keys())
#Rs_Tinit = np.zeros((len(factors), 3))
Rs_Tinit = np.zeros((len(factors), 2))
#rhos = np.zeros((len(factors), 3))
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

        Rs_Tinit[ii, 0] = max(Rs_Tinit[ii, 0], abs(R_PC1))
        Rs_Tinit[ii, 1] = max(Rs_Tinit[ii, 1], abs(R_PC2))
        #Rs_Tinit[ii, 2] = max(Rs_Tinit[ii, 2], abs(R_PC3))


Rs_Tinit = Rs_Tinit ** 2

#df_R_Tinit = pd.DataFrame(columns = ["PC1", "PC2", "PC3"], data=Rs_Tinit)
df_R_Tinit = pd.DataFrame(columns = ["PC1", "PC2"], data=Rs_Tinit)

df_R_Tinit['factor'] = factors
df_R_Tinit = df_R_Tinit.set_index("factor")
df_R_Tinit.to_csv("figure/CSV/explained_variance_df_R_Tinit.csv")
