# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 16:47:03 2021

@author: alex
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ks_2samp
import statsmodels.stats.multitest as multi
from scipy.stats import hypergeom


def sample_format(sample):
    new_name=[]
    lst=sample.columns.tolist()
    for k in range(len(lst)):
        split1=lst[k].split('r_')
        split2=split1[1].split('_tiret_')
        if len(split2)==1:
            counter=0
            for p in range (len(split2[0])-1):
                char1=split2[0][p]
                char2=split2[0][p+1]
                char=char1+char2
                if char=='_e' and p==len(split2[0])-2:
                    counter=1
                    split3=split2[0][:-2]
                    split3=split3+'[e]'
            if counter==0:
                split3=split2[0]           
        else:
            sting=split2[0]+split2[1]
            counter=0
            for p in range (len(sting)-1):
                char1=sting[p]
                char2=sting[p+1]
                char=char1+char2
                if char=='_e' and p==len(sting)-2:
                    counter=1
                    split3=sting[:-2]
                    split3=split3+'[e]'
            if counter==0:
                split3=sting
        new_name.append(split3)
    sample.columns = new_name

def model_dataframe(healthy_samples, Infected_samples):
    s_Healthy_rxns=healthy_samples.columns.tolist()
    s_INfected_rxns=Infected_samples.columns.tolist()
    INotH=[rxn for rxn in s_INfected_rxns if rxn not in s_Healthy_rxns]
    HNotI=[rxn for rxn in s_Healthy_rxns if rxn not in s_INfected_rxns]   
    for k in INotH:
        #print(k)
        healthy_samples.insert(healthy_samples.shape[1],k,0,False) 
    for k in HNotI:
        #print(k)
        Infected_samples.insert(Infected_samples.shape[1],k,0,False)
        
    data = pd.DataFrame(columns = ["reaction", "reactions_up_MW", "reactions_down_MW", "reactions_up_KS", "reactions_down_KS","p_ks","p_mw"]) 
    data_corrected = pd.DataFrame(columns = ["reaction", "reactions_up_MW", "reactions_down_MW", "reactions_up_KS", "reactions_down_KS","p_ks","p_mw"]) 
    
    p_mw=[]
    p_ks = []
    Fc=[]
    
    # p_value (Kolmogoro-smirnov et Mann-Whitney)
    s_Healthy_rxns=healthy_samples.columns.tolist()
    for k in s_Healthy_rxns:
        #print(k)
        rH=healthy_samples[[k]].values
        rH=np.ndarray.flatten(rH)
        rI=Infected_samples[[k]].values
        rI=np.ndarray.flatten(rI)
        t=[k for k in rH if k not in rI]
        if len(t)!=0:
            p_value=mannwhitneyu(rH,rI)[1]
        else:
            p_value=1#1
        p_mw.append(p_value)
        Fck=(np.mean(rI)-np.mean(rH))/(abs(np.mean(rI)+np.mean(rH)))
        Fc.append(Fck)
        p_value=ks_2samp(rH,rI)[1]
        p_ks.append(p_value)
    p_mw = np.array(p_mw)
    p_ks = np.array(p_ks)
    results_pw= multi.multipletests(p_mw, method = 'fdr_bh')[1]
    results_ks= multi.multipletests(p_ks, method = 'fdr_bh')[1]
    
    

    for k in range(len(p_ks)):
        data.loc[k]=[s_Healthy_rxns[k],0,0,0,0,p_ks[k],p_mw[k]]
        if p_ks[k]<0.05:
            if Fc[k]>=0.82:
                data.loc[k, ["reactions_up_KS"]] = 1
                data.loc[k, ["p_ks"]] = p_ks[k]              
            elif Fc[k]<=-0.82:
                data.loc[k, ["reactions_down_KS"]] = 1
                data.loc[k, ["p_ks"]] = p_ks[k]
                            
    for k in range(len(p_mw)):
        #df.loc[k]=[s_Healthy_rxns[k],0,0,0,0]
        if p_mw[k]<0.05:
            if Fc[k]>=0.82:
                data.loc[k, ["reactions_up_MW"]] = 1
                data.loc[k, ["p_mw"]] = p_mw[k]
            elif Fc[k]<=-0.82:
                data.loc[k, ["reactions_down_MW"]] = 1
                data.loc[k, ["p_mw"]] = p_mw[k]
                
    for k in range(len(results_ks)):
        data_corrected.loc[k]=[s_Healthy_rxns[k],0,0,0,0,results_ks[k],results_pw[k]]
        if results_ks[k]<0.05:
            if Fc[k]>=0.82:
                data_corrected.loc[k, ["reactions_up_KS"]] = 1
                data_corrected.loc[k, ["p_ks"]] = results_ks[k]
                
            elif Fc[k]<=-0.82:
                data_corrected.loc[k, ["reactions_down_KS"]] = 1
                data_corrected.loc[k, ["p_ks"]] = results_ks[k]
                
                
    for k in range(len(results_pw)):
        #df.loc[k]=[s_Healthy_rxns[k],0,0,0,0]
        if results_pw[k]<0.05:
            if Fc[k]>=0.82:
                data_corrected.loc[k, ["reactions_up_MW"]] = 1
                data_corrected.loc[k, ["p_mw"]] = results_pw[k]
            elif Fc[k]<=-0.82:
                data_corrected.loc[k, ["reactions_down_MW"]] = 1
                data_corrected.loc[k, ["p_mw"]] = results_pw[k]
    
    return data,data_corrected


def subsysteme_model(Healthy_model,Infected_model,data) :
    
    df_subsystems = pd.DataFrame(columns=["subsystem", "sizes", "reactions_up_MW", "reactions_down_MW", "reactions_up_KS", "reactions_down_KS"])  
    sub=[]
    number_reactions_up_MW = []
    number_reactions_down_MW =[]
    number_reactions_up_KS = []
    number_reactions_down_KS =[]
    number_reactions = []
    rxns_H=[]
    
    for k in range (len(Healthy_model.reactions)):
        rxns_H.append(Healthy_model.reactions[k]._id)
    
    rxns_I=[]
    for k in range (len(Infected_model.reactions)):
        rxns_I.append(Infected_model.reactions[k]._id)
    
    a=data['reaction'].to_list()
    for k in range (len(a)):
        if a[k] in rxns_H:
            index=rxns_H.index(a[k])
            sub_rxn=Healthy_model.reactions[index].subsystem
        if a[k] in rxns_I:
            index=rxns_I.index(a[k])
            sub_rxn=Infected_model.reactions[index].subsystem
        if sub_rxn not in sub :
            sub.append(sub_rxn)
            number_reactions.append(1)
            number_reactions_up_MW.append(data.reactions_up_MW.to_list()[k])
            number_reactions_down_MW.append(data.reactions_down_MW.to_list()[k])
            number_reactions_up_KS.append(data.reactions_up_KS.to_list()[k])
            number_reactions_down_KS.append(data.reactions_down_KS.to_list()[k])
            
        if sub_rxn in sub :
            idx=sub.index(sub_rxn)
            number_reactions[idx]=number_reactions[idx]+1
            number_reactions_up_MW[idx]+=data.reactions_up_MW.to_list()[k]
            number_reactions_down_MW[idx]+=data.reactions_down_MW.to_list()[k]
            number_reactions_up_KS[idx]+=data.reactions_up_KS.to_list()[k]
            number_reactions_down_KS[idx]+=data.reactions_down_KS.to_list()[k]
        
    for k in range (len(sub)):
        df_subsystems.loc[k]=[sub[k],number_reactions[k],number_reactions_up_MW[k],number_reactions_down_MW[k],number_reactions_up_KS[k],number_reactions_down_KS[k]]
    return df_subsystems
     
def enrichment_model(subsysteme_data):
    df_enrichment= pd.DataFrame(columns=["subsystem", "sizes", "up_MW", "down_MW", "up_KS", "down_KS"])  
    
    df_enrichment_corrected= pd.DataFrame(columns=["subsystem", "sizes", "up_MW", "down_MW", "up_KS", "down_KS"])  
    
    M=0
    for i in subsysteme_data.sizes.to_list():
        M+=i
    number_reactions_up_MW=0
    for i in subsysteme_data.reactions_up_MW.to_list():
        number_reactions_up_MW+=i
    
    number_reactions_down_MW=0
    for i in subsysteme_data.reactions_down_MW.to_list():
        number_reactions_down_MW+=i
        
    number_reactions_up_KS=0
    for i in subsysteme_data.reactions_up_KS.to_list():
        number_reactions_up_KS+=i
        
    number_reactions_down_KS=0
    for i in subsysteme_data.reactions_down_KS.to_list():
        number_reactions_down_KS+=i
        
    up_MW_cdf=[]
    down_MW_cdf=[]  
    up_KS_cdf=[]
    down_KS_cdf=[]  
    
        
    up_MW_sf=[]
    down_MW_sf=[]  
    up_KS_sf=[]
    down_KS_sf=[]  
    
    a=subsysteme_data.subsystem.to_list()

   
    for k in range (len(a)): 
        #M=subsystem_imat.reactions_up_MW[k] +  subsystem_imat.reactions_down_MW[k]
        #reaction Up MW
        N=subsysteme_data.sizes[k]
        n=number_reactions_up_MW
        
        #if n==0:
         #   x=0
        #else:
            #x=randint(1,n)
        
        x =subsysteme_data.reactions_up_MW[k]
        p=1-hypergeom.cdf(x,M, n, N)
        up_MW_cdf.append(p)
        u=hypergeom.sf(x,M, n, N)
        up_MW_sf.append(u)
        
        
        #reaction down MW
        N=subsysteme_data.sizes[k]
        n=number_reactions_down_MW
        x = subsysteme_data.reactions_down_MW[k] 
        p=1-hypergeom.cdf(x,M, n, N)
        down_MW_cdf.append(p)
        u=hypergeom.sf(x,M, n, N)
        down_MW_sf.append(u)
        
        #M=subsystem_imat.reactions_up_KS[k] +  subsystem_imat.reactions_down_KS[k]
        #reaction Up KS
        N=subsysteme_data.sizes[k]
        n=number_reactions_up_KS
        x =subsysteme_data.reactions_up_KS[k]
        p=1-hypergeom.cdf(x,M, n, N)
        up_KS_cdf.append(p)
        u=hypergeom.sf(x,M, n, N)
        up_KS_sf.append(u)
        
        
        #reaction down KS
        N=subsysteme_data.sizes[k]
        n=number_reactions_down_KS
        x =subsysteme_data.reactions_down_KS[k]
        p=1-hypergeom.cdf(x,M, n, N)
        down_KS_cdf.append(p)  
        u=hypergeom.sf(x,M, n, N)
        down_KS_sf.append(u)
    
    up_MW_cdf = np.array(up_MW_cdf)
    down_MW_cdf = np.array(down_MW_cdf)
    up_KS_cdf = np.array(up_KS_cdf)
    down_KS_cdf = np.array(down_KS_cdf)
    results_up_MW= multi.multipletests(up_MW_cdf, method = 'fdr_bh')[1]
    results_down_MW= multi.multipletests(down_MW_cdf, method = 'fdr_bh')[1]
    results_up_KS= multi.multipletests(up_KS_cdf, method = 'fdr_bh')[1]
    results_down_KS= multi.multipletests(down_KS_cdf, method = 'fdr_bh')[1]
    
    #for k in range (len(a)):
     #   df_enrichment.loc[k]=[subsysteme_data.subsystem[k],subsysteme_data.sizes[k],up_MW_cdf[k],down_MW_cdf[k],up_KS_cdf[k],down_KS_cdf[k]]
    
    for k in range (len(a)):
        df_enrichment_corrected.loc[k]=[subsysteme_data.subsystem[k],subsysteme_data.sizes[k],results_up_MW[k],results_down_MW[k],results_up_KS[k],results_down_KS[k]]
            
    return df_enrichment_corrected 

def altération(enrichment_corrected):
    df_altération= pd.DataFrame(columns=["subsystem", "sizes", "up_MW", "down_MW", "up_KS", "down_KS"])  
    for k in range(len(enrichment_corrected)):
         df_altération.loc[k]=[enrichment_corrected.subsystem[k],enrichment_corrected.sizes[k],0,0,0,0]
         if enrichment_corrected.up_MW[k]<0.05:
             df_altération.loc[k, ["up_MW"]] = 1
         if enrichment_corrected.down_MW[k]<0.05:
             df_altération.loc[k, ["down_MW"]] = 1
             
         if enrichment_corrected.up_KS[k]<0.05 and enrichment_corrected.up_KS[k]<enrichment_corrected.down_KS[k] :
             df_altération.loc[k, ["up_KS"]] = 1
         if enrichment_corrected.down_KS[k]<0.05 and enrichment_corrected.up_KS[k]>enrichment_corrected.down_KS[k] :
             df_altération.loc[k, ["down_KS"]] = 1
    return df_altération  

import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from os import listdir
myColors = ((53/255, 96/255, 149/255, 1.0), (1.0, 1.0, 1.0, 1.0), (167/255, 56/255, 44/255, 1.0))
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

def make_plot(df, encrichment_folder, file_name ):
    ax1 = sns.heatmap(df, cmap = cmap ,yticklabels=True, cbar=False,linewidths = 2)
    cbar = ax1.figure.colorbar(ax1.collections[0])
    cbar.set_ticks([-1, 0, 1])
    cbar.set_ticklabels(['down-regulated', 'unchanged', 'up-regulated'])
    ax1.xaxis.set_ticks_position('top')
    for item in ax1.get_yticklabels():
        item.set_rotation(30)
    for item in ax1.get_xticklabels():
        item.set_rotation(90)        
    fig = plt.gcf()
    if df.shape[1] > 2:
        fig.set_size_inches(2,(df.shape[0]/90)*20)
    else:
        fig.set_size_inches(1,(df.shape[0]/90)*20)
    
    plt.title(file_name)
    plt.savefig('data\\figure'+'\\'+file_name+'.png', bbox_inches = 'tight')
    plt.savefig('data\\figure'+'\\'+file_name+'.pdf', bbox_inches = 'tight')
    plt.show()
    plt.close()
 
def find_id(subsystem_list,subsystem):
    idx=[]
    for j in range(len(subsystem_list)):
        if subsystem_list[j]==subsystem:
            idx.append(j)
    return idx
    
    
    
def method_plot_MW(MEM):
    enrichment_folder = 'data\\subsystem_pvalue'
    columns =['data','subsystem','up_MW', 'down_MW','up_KS', 'down_KS']
    df = pd.DataFrame(columns=columns)
    
    files = listdir(enrichment_folder+"\\")
    
    for file in files:
            if "enrichment" in file and '.csv' in file and MEM in file:
                data = file.split("enrichment_")[1]
                data = data.split("_corrected_")[1]
                data = data.split(".csv")[0]
                df_tmp = pd.read_csv(enrichment_folder+"\\"+file, delimiter =";")
                df_tmp = df_tmp[['subsystem','up_MW', 'down_MW']]
                df_tmp['data'] = data
                df = df.append(df_tmp, ignore_index=True, sort=True)
                
    df['MW'] = [1 if u < 0.05 else -1 if d < 0.05 else 0 for u, d in list(zip(df.up_MW, df.down_MW))]      
    New=pd.DataFrame(columns=['subsystem']+list(set(df.data.tolist())))
    subsystem_list=df.subsystem.tolist()
    Data_list=df.data.tolist()
    MW_list=df.MW.tolist()
    subsystem_list_unique=list(set(df.subsystem.tolist()))
    cpt=0
    for k in subsystem_list_unique:
        idx=find_id(subsystem_list,k)
        value=[k,0,0,0,0,0]#order 293T,A549,CALU,Lung,NHBE
        for p in idx:
            dat=Data_list[p]
            if dat=='293T':
                value[1]=MW_list[p]
            if dat=='A549':
                  value[2]=MW_list[p]
            if dat=='CALU':
                 value[3]=MW_list[p]
            if dat=='Lung':
                 value[4]=MW_list[p]
            if dat=='NHBE':
                  value[5]=MW_list[p]
        c=0
        for t in range(1,len(value)):
            if value[t]!=0:
                c=1
        if c==1:        
            New.loc[cpt]=value
            cpt+=1
    string=MEM+' MW'
    New=New.set_index('subsystem')  
    New.to_csv(string+'.csv', header=True, index=True, sep=';', mode='w')
    Nez = pd.read_csv (string+'.csv',sep=';')
    Nez=Nez.set_index('subsystem')  
    make_plot(Nez, enrichment_folder,string)          
    
def method_plot_KS(MEM):
    enrichment_folder = 'data\\subsystem_pvalue' 


    columns =['data','subsystem','up_MW', 'down_MW','up_KS', 'down_KS']
    df = pd.DataFrame(columns=columns)
    
    files = listdir(enrichment_folder+"\\")
    
    for file in files:
            if "enrichment" in file and '.csv' in file and MEM in file:
                data = file.split("enrichment_")[1]
                data = data.split("_corrected_")[1]
                data = data.split(".csv")[0]
                df_tmp = pd.read_csv(enrichment_folder+"\\"+file, delimiter =";")
                df_tmp = df_tmp[['subsystem','up_KS', 'down_KS']]
                df_tmp['data'] = data
                df = df.append(df_tmp, ignore_index=True, sort=True)
    df['KS'] = [1 if u < 0.05 else -1 if d < 0.05 else 0 for u, d in list(zip(df.up_KS, df.down_KS))]      
    New=pd.DataFrame(columns=['subsystem']+list(set(df.data.tolist())))
    subsystem_list=df.subsystem.tolist()
    Data_list=df.data.tolist()
    KS_list=df.KS.tolist()
    subsystem_list_unique=list(set(df.subsystem.tolist()))
    cpt=0
    for k in subsystem_list_unique:

        idx=find_id(subsystem_list,k)

        value=[k,0,0,0,0,0]#order 293T,A549,CALU,Lung,NHBE
        for p in idx:
            dat=Data_list[p]
            if dat=='293T':
                value[1]=KS_list[p]
            if dat=='A549':
                  value[2]=KS_list[p]
            if dat=='CALU':
                 value[3]=KS_list[p]
            if dat=='Lung':
                 value[4]=KS_list[p]
            if dat=='NHBE':
                  value[5]=KS_list[p]
        c=0
        for t in range(1,len(value)):
            if value[t]!=0:
                c=1
        if c==1:        
            New.loc[cpt]=value
            cpt+=1
    string=MEM+' KS'
    New=New.set_index('subsystem')  
    New.to_csv(string+'.csv', header=True, index=True, sep=';', mode='w')
    Nez = pd.read_csv (string+'.csv',sep=';')
    Nez=Nez.set_index('subsystem')  


    make_plot(Nez, enrichment_folder,string)          
        
        
def MW_plot(data):



    enrichment_folder = 'data\\subsystem_pvalue' 
    columns = ['test','subsystem','up_MW', 'down_MW']
    df = pd.DataFrame(columns=columns)
    files = listdir(enrichment_folder+"\\")
    for file in files:
        if "enrichment" in file and '.csv' in file and data in file:
            method = file.split("enrichment_")[1]
            method = method.split("_corrected")[0]
            df_tmp = pd.read_csv(enrichment_folder+"\\"+file, sep=";")
            df_tmp = df_tmp[['subsystem','up_MW', 'down_MW']]
            df_tmp['method'] = method
            df = df.append(df_tmp, ignore_index=True, sort=True)    
    df['MW'] = [1 if u < 0.05 else -1 if d < 0.05 else 0 for u, d in list(zip(df.up_MW, df.down_MW))]
    New=pd.DataFrame(columns=['subsystem']+list(set(df.method.tolist())))
    subsystem_list=df.subsystem.tolist()
    Method_list=df.method.tolist()
    MW_list=df.MW.tolist()
    subsystem_list_unique=list(set(df.subsystem.tolist()))
    cpt=0
    for k in subsystem_list_unique:
        idx=find_id(subsystem_list,k)
        value=[k,0,0,0,0]#order Gimme, imat,init, Tinit
        for p in idx:
            dat=Method_list[p]
            if dat=='Gimme':
                value[1]=MW_list[p]
            if dat=='imat':
                  value[2]=MW_list[p]
            if dat=='init':
                 value[3]=MW_list[p]
            if dat=='Tinit':
                 value[4]=MW_list[p]
        c=0
        for t in range(1,len(value)):
            if value[t]!=0:
                c=1
        if c==1:        
            New.loc[cpt]=value
            cpt+=1
    string=data+' MW'
    New=New.set_index('subsystem')  
    New.to_csv(string+'.csv', header=True, index=True, sep=';', mode='w')
    Nez = pd.read_csv (string+'.csv',sep=';')
    Nez=Nez.set_index('subsystem')  


    make_plot(Nez, enrichment_folder,string)          
    
def KS_plot(data, min_common=0):
    enrichment_folder = 'data\\subsystem_pvalue' 
    columns = ['test','subsystem','up_KS', 'down_KS']
    df = pd.DataFrame(columns=columns)
    files = listdir(enrichment_folder+"\\")
    for file in files:
        if "enrichment" in file and '.csv' in file and data in file:
            method = file.split("enrichment_")[1]
            method = method.split("_corrected")[0]
            df_tmp = pd.read_csv(enrichment_folder+"\\"+file, sep=";")
            df_tmp = df_tmp[['subsystem','up_KS', 'down_KS']]
            df_tmp['method'] = method
            df = df.append(df_tmp, ignore_index=True, sort=True)   
    df['KS'] = [1 if u < 0.05 else -1 if d < 0.05 else 0 for u, d in list(zip(df.up_KS, df.down_KS))]
    New=pd.DataFrame(columns=['subsystem']+list(set(df.method.tolist())))
    subsystem_list=df.subsystem.tolist()
    Method_list=df.method.tolist()
    KS_list=df.KS.tolist()
    subsystem_list_unique=list(set(df.subsystem.tolist()))
    cpt=0
    for k in subsystem_list_unique:

        idx=find_id(subsystem_list,k)

        value=[k,0,0,0,0]#order Gimme, imat,init, Tinit
        for p in idx:
            dat=Method_list[p]
            if dat=='Gimme':
                value[1]=KS_list[p]
            if dat=='imat':
                 value[2]=KS_list[p]
            if dat=='init':
                 value[3]=KS_list[p]
            if dat=='Tinit':
                 value[4]=KS_list[p]

        c=0
        for t in range(1,len(value)):
            if value[t]!=0:
                c=1
        if c==1:        
            New.loc[cpt]=value
            cpt+=1
    string=data+' KS'
    New=New.set_index('subsystem')  
    New.to_csv(string+'.csv', header=True, index=True, sep=';', mode='w')
    Nez = pd.read_csv (string+'.csv',sep=';')
    Nez=Nez.set_index('subsystem')

    if min_common:
        Nez['enriched'] = np.sum(np.abs(Nez[Nez.columns[1:]].values), axis=0)
        Nez = Nez[Nez['enriched'] >= min_common]
        Nez = Nez.drop(columns=['enriched'])

    
    make_plot(Nez, enrichment_folder,string)  
 
