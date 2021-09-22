# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 13:32:07 2021

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

#PCA for Gimme
d=pd.read_csv('M_reaction_end_Gimme.csv',sep=',')
col_gimme=['Gimme_NHBE_H','Gimme_A549_H', 'Gimme_293T_H','Gimme_Lung_H','Gimme_CALU_H','Gimme_NHBE_I','Gimme_A549_I','Gimme_293T_I','Gimme_Lung_I','Gimme_CALU_I']
M_Gimme=d[col_gimme]
M_Gimme=np.transpose(M_Gimme.values)

data_Gimme_NHBE=['Gimme_NHBE_H','Gimme_NHBE_I']
data_Gimme_A549=['Gimme_A549_H','Gimme_A549_I']
data_Gimme_293T=['Gimme_293T_H','Gimme_293T_I']
data_Gimme_CALU=['Gimme_CALU_H','Gimme_CALU_I']
data_Gimme_Lung=['Gimme_Lung_H','Gimme_Lung_I']
data_Gimme=(data_Gimme_NHBE,data_Gimme_A549,data_Gimme_293T,data_Gimme_CALU,data_Gimme_Lung)

uninfected =col_gimme[:5]
infected=col_gimme[5:]
infection = (uninfected, infected)

groups = {"infection": infection, "Data":data_Gimme}
labels = {"infection": ("uninfected","infected"),
           "Data": ("NHBE","CALU","A549","293T","Lung")}


M_Gimme = M_Gimme[:,~np.all(M_Gimme==0, axis=0)]
M_Gimme= M_Gimme[:,~np.all(M_Gimme==1, axis=0)]
M_Gimme = np.array(M_Gimme, dtype=float)

M_Gimme = StandardScaler(with_mean=True, with_std=False).fit_transform(M_Gimme) #centering only

models=list(col_gimme)

#pca = PCA(n_components=3)
pca = PCA(n_components=2)
comps = pca.fit_transform(M_Gimme)
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
        plt.savefig("figure/Gimme_plot/Gimme_PC"+str(i1+1)+'_'+"PC"+str(i2+1)+'_'+group+".pdf",bbox_inches = 'tight')
        plt.show()
        
factors = list(groups.keys())
#Rs_gimme = np.zeros((len(factors), 3))
Rs_gimme = np.zeros((len(factors), 2))
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

        Rs_gimme[ii, 0] = max(Rs_gimme[ii, 0], abs(R_PC1))
        Rs_gimme[ii, 1] = max(Rs_gimme[ii, 1], abs(R_PC2))
        #Rs_gimme[ii, 2] = max(Rs_gimme[ii, 2], abs(R_PC3))


Rs_gimme = Rs_gimme ** 2

#df_R_gimme = pd.DataFrame(columns = ["PC1", "PC2", "PC3"], data=Rs_gimme)
df_R_gimme = pd.DataFrame(columns = ["PC1", "PC2"], data=Rs_gimme)

df_R_gimme['factor'] = factors
df_R_gimme = df_R_gimme.set_index("factor")
df_R_gimme.to_csv("figure/CSV/explained_variance_df_R_gimme.csv")
