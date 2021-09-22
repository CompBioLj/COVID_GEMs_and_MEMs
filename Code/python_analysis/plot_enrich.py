from pairwise_function import *

# Plot
data=['NHBE','A549','293T','Lung','CALU']
for k in data:
    KS_plot(k)
    MW_plot(k)
   

MEM=['Gimme','imat','Tinit','_init']
for k in MEM:
    method_plot_KS(k)
    method_plot_MW(k)
