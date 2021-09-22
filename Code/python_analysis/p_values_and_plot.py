# -*- coding: utf-8 -*-
"""
Created on Tue May 11 14:54:57 2021

@author: alex
"""

import cobra
import numpy as np
import pandas as pd
from cobra.sampling import sample
from os import listdir
from scipy.stats import mannwhitneyu, ks_2samp
import statsmodels.stats.multitest as multi
from pairwise_function import *
#from cobra.sampling import ACHRSampler

#293T
#Importing model
Tinit_293T_H=cobra.io.load_matlab_model('models/tINIT_model_H_293T_alt.mat')
Tinit_293T_I=cobra.io.load_matlab_model('models/tINIT_model_I_293T_alt.mat')

iMAT_293T_H=cobra.io.load_matlab_model('models/IMAT_model_H_293T_alt.mat')
iMAT_293T_I=cobra.io.load_matlab_model('models/IMAT_model_I_293T_alt.mat')

init_293T_H=cobra.io.load_matlab_model('models/INIT_model_H_293T_alt.mat')
init_293T_I=cobra.io.load_matlab_model('models/INIT_model_I_293T_alt.mat')

Gimme_293T_H=cobra.io.load_matlab_model('models/Gimme_model_H_293T_alt.mat')
Gimme_293T_I=cobra.io.load_matlab_model('models/Gimme_model_I_293T_alt.mat')

#Sampling 

s_Healthy_Tinit_293T = pd.read_csv ('flux_sampling/Tinit_293T_sample_H.csv')
s_Infected_Tinit_293T = pd.read_csv ('flux_sampling/Tinit_293T_sample_I.csv')
s_Healthy_iMAT_293T = pd.read_csv ('flux_sampling/iMAT_293T_sample_H.csv')
s_Infected_iMAT_293T = pd.read_csv ('flux_sampling/iMAT_293T_sample_I.csv')
s_Healthy_init_293T=pd.read_csv ('flux_sampling/init_293T_sample_H.csv')
s_Infected_init_293T=pd.read_csv ('flux_sampling/init_293T_sample_I.csv')
s_Healthy_Gimme_293T=pd.read_csv ('flux_sampling/Gimme_293T_sample_H.csv')
s_Infected_Gimme_293T=pd.read_csv ('flux_sampling/Gimme_293T_sample_I.csv')

#_A549
#Importing model
Tinit_A549_H=cobra.io.load_matlab_model('models/tINIT_model_H_A549_alt.mat')
Tinit_A549_I=cobra.io.load_matlab_model('models/tINIT_model_I_A549_alt.mat')

iMAT_A549_H=cobra.io.load_matlab_model('models/IMAT_model_H_A549_alt.mat')
iMAT_A549_I=cobra.io.load_matlab_model('models/iMAT_model_I_A549_alt.mat')

init_A549_H=cobra.io.load_matlab_model('models/INIT_model_H_A549_alt.mat')
init_A549_I=cobra.io.load_matlab_model('models/INIT_model_I_A549_alt.mat')

Gimme_A549_H=cobra.io.load_matlab_model('models/Gimme_model_H_A549_alt.mat')
Gimme_A549_I=cobra.io.load_matlab_model('models/Gimme_model_I_A549_alt.mat')

#Sampling 
s_Healthy_Tinit_A549 = pd.read_csv ('flux_sampling/Tinit_A549_sample_H.csv')
s_Infected_Tinit_A549 = pd.read_csv ('flux_sampling/Tinit_A549_sample_I.csv')
s_Healthy_iMAT_A549 = pd.read_csv ('flux_sampling/iMAT_A549_sample_H.csv')
s_Infected_iMAT_A549 = pd.read_csv ('flux_sampling/iMAT_A549_sample_I.csv')
s_Healthy_init_A549=pd.read_csv ('flux_sampling/init_A549_sample_H.csv')
s_Infected_init_A549=pd.read_csv ('flux_sampling/init_A549_sample_I.csv')
s_Healthy_Gimme_A549=pd.read_csv ('flux_sampling/Gimme_A549_sample_H.csv')
s_Infected_Gimme_A549=pd.read_csv ('flux_sampling/Gimme_A549_sample_I.csv')
#_CALU
#Importing model
Tinit_CALU_H=cobra.io.load_matlab_model('models/tINIT_model_H_CALU_alt.mat')
Tinit_CALU_I=cobra.io.load_matlab_model('models/tINIT_model_I_CALU_alt.mat')
iMAT_CALU_H=cobra.io.load_matlab_model('models/IMAT_model_H_CALU_alt.mat')
iMAT_CALU_I=cobra.io.load_matlab_model('models/iMAT_model_I_CALU_alt.mat')
init_CALU_H=cobra.io.load_matlab_model('models/INIT_model_H_CALU_alt.mat')
init_CALU_I=cobra.io.load_matlab_model('models/INIT_model_I_CALU_alt.mat')
Gimme_CALU_H=cobra.io.load_matlab_model('models/Gimme_model_H_CALU_alt.mat')
Gimme_CALU_I=cobra.io.load_matlab_model('models/Gimme_model_I_CALU_alt.mat')
#Sampling 
s_Healthy_Tinit_CALU = pd.read_csv ('flux_sampling/Tinit_CALU_sample_H.csv')
s_Infected_Tinit_CALU = pd.read_csv ('flux_sampling/Tinit_CALU_sample_I.csv')
s_Healthy_iMAT_CALU = pd.read_csv ('flux_sampling/iMAT_CALU_sample_H.csv')
s_Infected_iMAT_CALU = pd.read_csv ('flux_sampling/iMAT_CALU_sample_I.csv')
s_Healthy_init_CALU=pd.read_csv ('flux_sampling/init_CALU_sample_H.csv')
s_Infected_init_CALU=pd.read_csv ('flux_sampling/init_CALU_sample_I.csv')
s_Healthy_Gimme_CALU=pd.read_csv ('flux_sampling/Gimme_CALU_sample_H.csv')
s_Infected_Gimme_CALU=pd.read_csv ('flux_sampling/Gimme_CALU_sample_I.csv')
#_Lung
#Importing model
Tinit_Lung_H=cobra.io.load_matlab_model('models/tINIT_model_H_Lung_alt.mat')
Tinit_Lung_I=cobra.io.load_matlab_model('models/tINIT_model_I_Lung_alt.mat')
iMAT_Lung_H=cobra.io.load_matlab_model('models/IMAT_model_H_Lung_alt.mat')
iMAT_Lung_I=cobra.io.load_matlab_model('models/iMAT_model_I_Lung_alt.mat')
init_Lung_H=cobra.io.load_matlab_model('models/INIT_model_H_Lung_alt.mat')
init_Lung_I=cobra.io.load_matlab_model('models/INIT_model_I_Lung_alt.mat')
Gimme_Lung_H=cobra.io.load_matlab_model('models/Gimme_model_H_Lung_alt.mat')
Gimme_Lung_I=cobra.io.load_matlab_model('models/Gimme_model_I_Lung_alt.mat')
#Sampling 
s_Healthy_Tinit_Lung = pd.read_csv ('flux_sampling/Tinit_Lung_sample_H.csv')
s_Infected_Tinit_Lung = pd.read_csv ('flux_sampling/Tinit_Lung_sample_I.csv')
s_Healthy_iMAT_Lung = pd.read_csv ('flux_sampling/iMAT_Lung_sample_H.csv')
s_Infected_iMAT_Lung = pd.read_csv ('flux_sampling/iMAT_Lung_sample_I.csv')
s_Healthy_init_Lung=pd.read_csv ('flux_sampling/init_Lung_sample_H.csv')
s_Infected_init_Lung=pd.read_csv ('flux_sampling/init_Lung_sample_I.csv')
s_Healthy_Gimme_Lung=pd.read_csv ('flux_sampling/Gimme_Lung_sample_H.csv')
s_Infected_Gimme_Lung=pd.read_csv ('flux_sampling/Gimme_Lung_sample_I.csv')
#_NHBE
#Importing model
Tinit_NHBE_H=cobra.io.load_matlab_model('models/tINIT_model_H_NHBE.mat')
Tinit_NHBE_I=cobra.io.load_matlab_model('models/tINIT_model_I_NHBE.mat')
iMAT_NHBE_H=cobra.io.load_matlab_model('models/IMAT_model_H_NHBE_alt.mat')
iMAT_NHBE_I=cobra.io.load_matlab_model('models/iMAT_model_I_NHBE_alt.mat')
init_NHBE_H=cobra.io.load_matlab_model('models/INIT_model_H_NHBE.mat')
init_NHBE_I=cobra.io.load_matlab_model('models/INIT_model_I_NHBE.mat')
Gimme_NHBE_H=cobra.io.load_matlab_model('models/Gimme_model_H_NHBE_alt.mat')
Gimme_NHBE_I=cobra.io.load_matlab_model('models/Gimme_model_I_NHBE_alt.mat')
#Sampling 
s_Healthy_Tinit_NHBE = pd.read_csv ('flux_sampling/Tinit_NHBE_sample_H.csv')
s_Infected_Tinit_NHBE = pd.read_csv ('flux_sampling/Tinit_NHBE_sample_I.csv')
s_Healthy_iMAT_NHBE = pd.read_csv ('flux_sampling/iMAT_NHBE_sample_H.csv')
s_Infected_iMAT_NHBE = pd.read_csv ('flux_sampling/iMAT_NHBE_sample_I.csv')
s_Healthy_init_NHBE=pd.read_csv ('flux_sampling/init_NHBE_sample_H.csv')
s_Infected_init_NHBE=pd.read_csv ('flux_sampling/init_NHBE_sample_I.csv')
s_Healthy_Gimme_NHBE=pd.read_csv ('flux_sampling/Gimme_NHBE_sample_H.csv')
s_Infected_Gimme_NHBE=pd.read_csv ('flux_sampling/Gimme_NHBE_sample_I.csv')

 
    
sample_format(s_Healthy_Tinit_293T)
sample_format(s_Infected_Tinit_293T)
sample_format(s_Healthy_iMAT_293T) 
sample_format(s_Infected_iMAT_293T)
sample_format(s_Healthy_init_293T)
sample_format(s_Infected_init_293T)
sample_format(s_Healthy_Gimme_293T)
sample_format(s_Infected_Gimme_293T)

sample_format(s_Healthy_Tinit_A549)
sample_format(s_Infected_Tinit_A549)
sample_format(s_Healthy_iMAT_A549) 
sample_format(s_Infected_iMAT_A549)
sample_format(s_Healthy_init_A549)
sample_format(s_Infected_init_A549)
sample_format(s_Healthy_Gimme_A549)
sample_format(s_Infected_Gimme_A549)
 
sample_format(s_Healthy_Tinit_CALU)
sample_format(s_Infected_Tinit_CALU)
sample_format(s_Healthy_iMAT_CALU) 
sample_format(s_Infected_iMAT_CALU)
sample_format(s_Healthy_init_CALU)
sample_format(s_Infected_init_CALU)
sample_format(s_Healthy_Gimme_CALU)
sample_format(s_Infected_Gimme_CALU)

sample_format(s_Healthy_Tinit_Lung)
sample_format(s_Infected_Tinit_Lung)
sample_format(s_Healthy_iMAT_Lung) 
sample_format(s_Infected_iMAT_Lung)
sample_format(s_Healthy_init_Lung)
sample_format(s_Infected_init_Lung)
sample_format(s_Healthy_Gimme_Lung)
sample_format(s_Infected_Gimme_Lung)
 
sample_format(s_Healthy_Tinit_NHBE)
sample_format(s_Infected_Tinit_NHBE)
sample_format(s_Healthy_iMAT_NHBE) 
sample_format(s_Infected_iMAT_NHBE)
sample_format(s_Healthy_init_NHBE)
sample_format(s_Infected_init_NHBE)
sample_format(s_Healthy_Gimme_NHBE)
sample_format(s_Infected_Gimme_NHBE)

Tinit_data_293T,Tinit_data_corrected_293T=model_dataframe(s_Healthy_Tinit_293T, s_Infected_Tinit_293T)
iMAT_data_293T,iMAT_data_corrected_293T=model_dataframe(s_Healthy_iMAT_293T,s_Infected_iMAT_293T)
init_data_293T,init_data_corrected_293T=model_dataframe(s_Healthy_init_293T,s_Infected_init_293T)
Gimme_data_293T,Gimme_data_corrected_293T=model_dataframe(s_Healthy_Gimme_293T,s_Infected_Gimme_293T)

# Save in txt file corrected p values 
Tinit_data_corrected_293T.to_csv('data/reaction_pvalue/Tinit_data_corrected_293T.csv', header=True, index=None, sep=';', mode='w')

iMAT_data_corrected_293T.to_csv('data/reaction_pvalue/iMAT_data_corrected_293T.csv', header=True, index=None, sep=';', mode='w')

init_data_corrected_293T.to_csv('data/reaction_pvalue/init_data_corrected_293T.csv', header=True, index=None, sep=';', mode='w')

Gimme_data_corrected_293T.to_csv('data/reaction_pvalue/Gimme_data_corrected_293T.csv', header=True, index=None, sep=';', mode='w')
#A549
Tinit_data_A549,Tinit_data_corrected_A549=model_dataframe(s_Healthy_Tinit_A549, s_Infected_Tinit_A549)
iMAT_data_A549,iMAT_data_corrected_A549=model_dataframe(s_Healthy_iMAT_A549,s_Infected_iMAT_A549)
init_data_A549,init_data_corrected_A549=model_dataframe(s_Healthy_init_A549,s_Infected_init_A549)
Gimme_data_A549,Gimme_data_corrected_A549=model_dataframe(s_Healthy_Gimme_A549,s_Infected_Gimme_A549)

# Save in txt file corrected p values 
Tinit_data_corrected_A549.to_csv('data/reaction_pvalue/Tinit_data_corrected_A549.csv', header=True, index=None, sep=';', mode='w')

iMAT_data_corrected_A549.to_csv('data/reaction_pvalue/iMAT_data_corrected_A549.csv', header=True, index=None, sep=';', mode='w')

init_data_corrected_A549.to_csv('data/reaction_pvalue/init_data_corrected_A549.csv', header=True, index=None, sep=';', mode='w')

Gimme_data_corrected_A549.to_csv('data/reaction_pvalue/Gimme_data_corrected_A549.csv', header=True, index=None, sep=';', mode='w')

#CALU
Tinit_data_CALU,Tinit_data_corrected_CALU=model_dataframe(s_Healthy_Tinit_CALU, s_Infected_Tinit_CALU)
iMAT_data_CALU,iMAT_data_corrected_CALU=model_dataframe(s_Healthy_iMAT_CALU,s_Infected_iMAT_CALU)
init_data_CALU,init_data_corrected_CALU=model_dataframe(s_Healthy_init_CALU,s_Infected_init_CALU)
Gimme_data_CALU,Gimme_data_corrected_CALU=model_dataframe(s_Healthy_Gimme_CALU,s_Infected_Gimme_CALU)

# Save in txt file corrected p values 
Tinit_data_corrected_CALU.to_csv('data/reaction_pvalue/Tinit_data_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')

iMAT_data_corrected_CALU.to_csv('data/reaction_pvalue/iMAT_data_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')

init_data_corrected_CALU.to_csv('data/reaction_pvalue/init_data_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')

Gimme_data_corrected_CALU.to_csv('data/reaction_pvalue/Gimme_data_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')

#Lung
Tinit_data_Lung,Tinit_data_corrected_Lung=model_dataframe(s_Healthy_Tinit_Lung, s_Infected_Tinit_Lung)
iMAT_data_Lung,iMAT_data_corrected_Lung=model_dataframe(s_Healthy_iMAT_Lung,s_Infected_iMAT_Lung)
init_data_Lung,init_data_corrected_Lung=model_dataframe(s_Healthy_init_Lung,s_Infected_init_Lung)
Gimme_data_Lung,Gimme_data_corrected_Lung=model_dataframe(s_Healthy_Gimme_Lung,s_Infected_Gimme_Lung)

# Save in txt file corrected p values 
Tinit_data_corrected_Lung.to_csv('data/reaction_pvalue/Tinit_data_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')

iMAT_data_corrected_Lung.to_csv('data/reaction_pvalue/iMAT_data_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')

init_data_corrected_Lung.to_csv('data/reaction_pvalue/init_data_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')

Gimme_data_corrected_Lung.to_csv('data/reaction_pvalue/Gimme_data_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')
 

 #NHBE
Tinit_data_NHBE,Tinit_data_corrected_NHBE=model_dataframe(s_Healthy_Tinit_NHBE, s_Infected_Tinit_NHBE)
iMAT_data_NHBE,iMAT_data_corrected_NHBE=model_dataframe(s_Healthy_iMAT_NHBE,s_Infected_iMAT_NHBE)
init_data_NHBE,init_data_corrected_NHBE=model_dataframe(s_Healthy_init_NHBE,s_Infected_init_NHBE)
Gimme_data_NHBE,Gimme_data_corrected_NHBE=model_dataframe(s_Healthy_Gimme_NHBE,s_Infected_Gimme_NHBE)

# Save in txt file corrected p values 
Tinit_data_corrected_NHBE.to_csv('data/reaction_pvalue/Tinit_data_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')

iMAT_data_corrected_NHBE.to_csv('data/reaction_pvalue/iMAT_data_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')

init_data_corrected_NHBE.to_csv('data/reaction_pvalue/init_data_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')

Gimme_data_corrected_NHBE.to_csv('data/reaction_pvalue/Gimme_data_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')

#293T
subsystem_imat_corrected_293T=subsysteme_model(iMAT_293T_H,iMAT_293T_I,iMAT_data_corrected_293T) 
subsystem_Tinit_corrected_293T=subsysteme_model(Tinit_293T_H,Tinit_293T_I,Tinit_data_corrected_293T) 
subsystem_init_corrected_293T=subsysteme_model(init_293T_H,init_293T_I,init_data_corrected_293T)  
subsystem_Gimme_corrected_293T=subsysteme_model(Gimme_293T_H,Gimme_293T_I,Gimme_data_corrected_293T)
enrichment_Gimme_corrected_293T=enrichment_model(subsystem_Gimme_corrected_293T)

enrichment_Gimme_corrected_293T.to_csv('data/subsystem_pvalue/enrichment_Gimme_corrected_293T.csv', header=True, index=None, sep=';', mode='w')
enrichment_imat_corrected_293T=enrichment_model(subsystem_imat_corrected_293T) 
enrichment_Tinit_corrected_293T=enrichment_model(subsystem_Tinit_corrected_293T) 
enrichment_init_corrected_293T=enrichment_model(subsystem_init_corrected_293T) 
enrichment_imat_corrected_293T.to_csv('data/subsystem_pvalue/enrichment_imat_corrected_293T.csv', header=True, index=None, sep=';', mode='w')
enrichment_Tinit_corrected_293T.to_csv('data/subsystem_pvalue/enrichment_Tinit_corrected_293T.csv', header=True, index=None, sep=';', mode='w')
enrichment_init_corrected_293T.to_csv('data/subsystem_pvalue/enrichment_init_corrected_293T.csv', header=True, index=None, sep=';', mode='w')

#_A549
subsystem_Gimme_corrected_A549=subsysteme_model(Gimme_A549_H,Gimme_A549_I,Gimme_data_corrected_A549)
enrichment_Gimme_corrected_A549=enrichment_model(subsystem_Gimme_corrected_A549)

enrichment_Gimme_corrected_A549.to_csv('data/subsystem_pvalue/enrichment_Gimme_corrected_A549.csv', header=True, index=None, sep=';', mode='w')
subsystem_imat_corrected_A549=subsysteme_model(iMAT_A549_H,iMAT_A549_I,iMAT_data_corrected_A549) 
subsystem_Tinit_corrected_A549=subsysteme_model(Tinit_A549_H,Tinit_A549_I,Tinit_data_corrected_A549) 
subsystem_init_corrected_A549=subsysteme_model(init_A549_H,init_A549_I,init_data_corrected_A549) 

enrichment_imat_corrected_A549=enrichment_model(subsystem_imat_corrected_A549) 
enrichment_Tinit_corrected_A549=enrichment_model(subsystem_Tinit_corrected_A549) 
enrichment_init_corrected_A549=enrichment_model(subsystem_init_corrected_A549) 

enrichment_imat_corrected_A549.to_csv('data/subsystem_pvalue/enrichment_imat_corrected_A549.csv', header=True, index=None, sep=';', mode='w')
enrichment_Tinit_corrected_A549.to_csv('data/subsystem_pvalue/enrichment_Tinit_corrected_A549.csv', header=True, index=None, sep=';', mode='w')
enrichment_init_corrected_A549.to_csv('data/subsystem_pvalue/enrichment_init_corrected_A549.csv', header=True, index=None, sep=';', mode='w') 

#_Lung
subsystem_Gimme_corrected_Lung=subsysteme_model(Gimme_Lung_H,Gimme_Lung_I,Gimme_data_corrected_Lung)
enrichment_Gimme_corrected_Lung=enrichment_model(subsystem_Gimme_corrected_Lung)

enrichment_Gimme_corrected_Lung.to_csv('data/subsystem_pvalue/enrichment_Gimme_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')

subsystem_imat_corrected_Lung=subsysteme_model(iMAT_Lung_H,iMAT_Lung_I,iMAT_data_corrected_Lung) 
subsystem_Tinit_corrected_Lung=subsysteme_model(Tinit_Lung_H,Tinit_Lung_I,Tinit_data_corrected_Lung) 
subsystem_init_corrected_Lung=subsysteme_model(init_Lung_H,init_Lung_I,init_data_corrected_Lung) 
enrichment_imat_corrected_Lung=enrichment_model(subsystem_imat_corrected_Lung) 
enrichment_Tinit_corrected_Lung=enrichment_model(subsystem_Tinit_corrected_Lung) 
enrichment_init_corrected_Lung=enrichment_model(subsystem_init_corrected_Lung) 
enrichment_imat_corrected_Lung.to_csv('data/subsystem_pvalue/enrichment_imat_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')
enrichment_Tinit_corrected_Lung.to_csv('data/subsystem_pvalue/enrichment_Tinit_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')
enrichment_init_corrected_Lung.to_csv('data/subsystem_pvalue/enrichment_init_corrected_Lung.csv', header=True, index=None, sep=';', mode='w') 
 
 #_NHBE
subsystem_Gimme_corrected_NHBE=subsysteme_model(Gimme_NHBE_H,Gimme_NHBE_I,Gimme_data_corrected_NHBE)
enrichment_Gimme_corrected_NHBE=enrichment_model(subsystem_Gimme_corrected_NHBE)

enrichment_Gimme_corrected_NHBE.to_csv('data/subsystem_pvalue/enrichment_Gimme_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')
subsystem_imat_corrected_NHBE=subsysteme_model(iMAT_NHBE_H,iMAT_NHBE_I,iMAT_data_corrected_NHBE) 
subsystem_Tinit_corrected_NHBE=subsysteme_model(Tinit_NHBE_H,Tinit_NHBE_I,Tinit_data_corrected_NHBE) 
subsystem_init_corrected_NHBE=subsysteme_model(init_NHBE_H,init_NHBE_I,init_data_corrected_NHBE) 
enrichment_imat_corrected_NHBE=enrichment_model(subsystem_imat_corrected_NHBE) 
enrichment_Tinit_corrected_NHBE=enrichment_model(subsystem_Tinit_corrected_NHBE) 
enrichment_init_corrected_NHBE=enrichment_model(subsystem_init_corrected_NHBE) 
enrichment_imat_corrected_NHBE.to_csv('data/subsystem_pvalue/enrichment_imat_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')
enrichment_Tinit_corrected_NHBE.to_csv('data/subsystem_pvalue/enrichment_Tinit_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')
enrichment_init_corrected_NHBE.to_csv('data/subsystem_pvalue/enrichment_init_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')
 #_CALU
subsystem_Gimme_corrected_CALU=subsysteme_model(Gimme_CALU_H,Gimme_CALU_I,Gimme_data_corrected_CALU)
enrichment_Gimme_corrected_CALU=enrichment_model(subsystem_Gimme_corrected_CALU)

enrichment_Gimme_corrected_CALU.to_csv('data/subsystem_pvalue/enrichment_Gimme_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')
subsystem_imat_corrected_CALU=subsysteme_model(iMAT_CALU_H,iMAT_CALU_I,iMAT_data_corrected_CALU) 
subsystem_Tinit_corrected_CALU=subsysteme_model(Tinit_CALU_H,Tinit_CALU_I,Tinit_data_corrected_CALU) 
subsystem_init_corrected_CALU=subsysteme_model(init_CALU_H,init_CALU_I,init_data_corrected_CALU) 

enrichment_imat_corrected_CALU=enrichment_model(subsystem_imat_corrected_CALU) 
enrichment_Tinit_corrected_CALU=enrichment_model(subsystem_Tinit_corrected_CALU) 
enrichment_init_corrected_CALU=enrichment_model(subsystem_init_corrected_CALU) 
enrichment_imat_corrected_CALU.to_csv('data/subsystem_pvalue/enrichment_imat_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')
enrichment_Tinit_corrected_CALU.to_csv('data/subsystem_pvalue/enrichment_Tinit_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')
enrichment_init_corrected_CALU.to_csv('data/subsystem_pvalue/enrichment_init_corrected_CALU.csv', header=True, index=None, sep=';', mode='w') 
#293T            
altération_Gimme_corrected_293T=altération(enrichment_Gimme_corrected_293T) 
altération_Gimme_corrected_293T.to_csv('data/subsystem_binary/altération_Gimme_corrected_293T.csv', header=True, index=None, sep=';', mode='w')
altération_imat_corrected_293T=altération(enrichment_imat_corrected_293T) 
altération_Tinit_corrected_293T=altération(enrichment_Tinit_corrected_293T) 
altération_init_corrected_293T=altération(enrichment_init_corrected_293T) 
altération_imat_corrected_293T.to_csv('data/subsystem_binary/altération_imat_corrected_293T.csv', header=True, index=None, sep=';', mode='w')
altération_Tinit_corrected_293T.to_csv('data/subsystem_binary/altération_Tinit_corrected_293T.csv', header=True, index=None, sep=';', mode='w')
altération_init_corrected_293T.to_csv('data/subsystem_binary/altération_init_corrected_293T.csv', header=True, index=None, sep=';', mode='w')

#_A549     
altération_Gimme_corrected_A549=altération(enrichment_Gimme_corrected_A549) 
altération_Gimme_corrected_A549.to_csv('data/subsystem_binary/altération_Gimme_corrected_A549.csv', header=True, index=None, sep=';', mode='w')       
altération_imat_corrected_A549=altération(enrichment_imat_corrected_A549) 
altération_Tinit_corrected_A549=altération(enrichment_Tinit_corrected_A549) 
altération_init_corrected_A549=altération(enrichment_init_corrected_A549) 
altération_imat_corrected_A549.to_csv('data/subsystem_binary/altération_imat_corrected_A549.csv', header=True, index=None, sep=';', mode='w')
altération_Tinit_corrected_A549.to_csv('data/subsystem_binary/altération_Tinit_corrected_A549.csv', header=True, index=None, sep=';', mode='w')
altération_init_corrected_A549.to_csv('data/subsystem_binary/altération_init_corrected_A549.csv', header=True, index=None, sep=';', mode='w')
#_Lung  
altération_Gimme_corrected_Lung=altération(enrichment_Gimme_corrected_Lung) 
altération_Gimme_corrected_Lung.to_csv('data/subsystem_binary/altération_Gimme_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')          
altération_imat_corrected_Lung=altération(enrichment_imat_corrected_Lung) 
altération_Tinit_corrected_Lung=altération(enrichment_Tinit_corrected_Lung) 
altération_init_corrected_Lung=altération(enrichment_init_corrected_Lung) 
altération_imat_corrected_Lung.to_csv('data/subsystem_binary/altération_imat_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')
altération_Tinit_corrected_Lung.to_csv('data/subsystem_binary/altération_Tinit_corrected_Lung.csv', header=True, index=None, sep=';', mode='w')
altération_init_corrected_Lung.to_csv('data/subsystem_binary/altération_init_corrected_Lung.csv', header=True, index=None, sep=';', mode='w') 
#_NHBE   
altération_Gimme_corrected_NHBE=altération(enrichment_Gimme_corrected_NHBE) 
altération_Gimme_corrected_NHBE.to_csv('data/subsystem_binary/altération_Gimme_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')         
altération_imat_corrected_NHBE=altération(enrichment_imat_corrected_NHBE) 
altération_Tinit_corrected_NHBE=altération(enrichment_Tinit_corrected_NHBE) 
altération_init_corrected_NHBE=altération(enrichment_init_corrected_NHBE) 
altération_imat_corrected_NHBE.to_csv('data/subsystem_binary/altération_imat_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')
altération_Tinit_corrected_NHBE.to_csv('data/subsystem_binary/altération_Tinit_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w')
altération_init_corrected_NHBE.to_csv('data/subsystem_binary/altération_init_corrected_NHBE.csv', header=True, index=None, sep=';', mode='w') 
 #_CALU      
altération_Gimme_corrected_CALU=altération(enrichment_Gimme_corrected_CALU) 
altération_Gimme_corrected_CALU.to_csv('data/subsystem_binary/altération_Gimme_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')   
altération_imat_corrected_CALU=altération(enrichment_imat_corrected_CALU) 
altération_Tinit_corrected_CALU=altération(enrichment_Tinit_corrected_CALU) 
altération_init_corrected_CALU=altération(enrichment_init_corrected_CALU) 
altération_imat_corrected_CALU.to_csv('data/subsystem_binary/altération_imat_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')
altération_Tinit_corrected_CALU.to_csv('data/subsystem_binary/altération_Tinit_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')
altération_init_corrected_CALU.to_csv('data/subsystem_binary/altération_init_corrected_CALU.csv', header=True, index=None, sep=';', mode='w')

# Plot

data=['NHBE','A549','293T','Lung','CALU']
for k in data:
    KS_plot(k)
    MW_plot(k)
   

MEM=['Gimme','imat','Tinit','_init']
for k in MEM:
    method_plot_KS(k)
    method_plot_MW(k)
