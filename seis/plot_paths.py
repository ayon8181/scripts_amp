import imp
import math
import numpy as np
import matplotlib.pyplot as plt
import SeisTomoPy

model = 'S40RTS'  
# Parameter to be plotted
para = 'VS'    
# Latitutde of the starting point of the cross section
elat = -60    
# Longitude of the starting point of the cross section
elon = -49   
# Latitutde of the ending point of the cross section
slat = 20 
# Longitude of the ending point of the cross section
slon = 69  
# Maximal velocity perturbations for the colorbar
Vmax = 2             
# List of seismic phases
phlist = 'SSS Sdiff'
# 1-D model to compute the paths
# Can be one of the default models 
model1d = 'prem'
# Position of stations and event
EVT = np.array([10,100])#np.loadtxt('/home1/09038/ayon8181/SeisTomoPy_V3/Documentation/SeisTomoPy_notebook/files/event.xy')
STA = np.array([120,0])#lsnp.loadtxt('/home1/09038/ayon8181//SeisTomoPy_V3/Documentation/SeisTomoPy_notebook/files/station.xy')

# Running path_plot
fig=SeisTomoPy.path_plot(model,para,Vmax,elat,elon,slat,slon,EVT,STA,phlist,model1d)
plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/paths2.png")