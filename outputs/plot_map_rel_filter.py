#!/usr/bin/env python
import numpy as np
import csv
import pygmt
import pandas as pd
import os
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
import math
import multiprocessing as mp
import csv
#import seaborn as sns
import scipy.stats as stats
import geopy.distance
from geographiclib.geodesic import Geodesic
from matplotlib.lines import Line2D
#from mpi4py import MPI

file_names      = ["_real","_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['darkred','cyan','yellowgreen','chocolate','black']    
lss   =[2,2.5,2,1.5,1] 
style =['None','solid',"solid","solid","solid"]
labels=["real_data",'S40RTS_1D_Crust','S40RTS_3D_crust',"GLAD_M25","PREM_Crust2.0"] #"PREM_QRFSI12",
fccs  =['hotpink','None','None','None',"None"] 
markers=[".","+","o","*","x"]
alp   =[0.5,1,1,1,1]
df=[]

plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/plots_filter_rel"

ev_list=[]
with open("./../event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])

ph="SS"
k=1
moho={}
with open("/scratch1/09038/ayon8181/scripts_amp/outputs/depthtomoho.xyz",'r') as txt:
    data = csv.reader(txt, skipinitialspace=True, delimiter=" ")
    for row in data:
        lats   = math.floor(float(row[1]))
        lons   = math.floor(float(row[0]))
        if lats not in moho.keys():
            moho[lats] = {}
        moho[lats][lons] = float(row[2])

def points(evla,evlo,evdp,stla,stlo,gcarc,lim1,lim2):
    
    az=Geodesic.WGS84.Inverse(evla, evlo, stla, stlo)['azi1']
    st_di=gcarc*lim1
    start=geopy.distance.distance(nautical=60*st_di).destination((evla,evlo),bearing=az)
    m_di=gcarc*(1/2.0)
    mid=geopy.distance.distance(nautical=60*m_di).destination((evla,evlo),bearing=az)
    en_di=gcarc*lim2
    end=geopy.distance.distance(nautical=60*en_di).destination((evla,evlo),bearing=az)
    #print(start[1],float("{:3.2f}".format(start.latitude)),az)
    
    return [[start[1],start[0]],[mid[1],mid[0]],[end[1],end[0]]]

df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_plot=df[df[str(16+k*7+5)+"_real"].notna()]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]
#df_plot=df_plot[(np.abs(df_plot[str(16+k*7+2)+"_real"])<15) & (np.abs(df_plot[str(16+k*7+2)+"_1D"])<15) &(np.abs(df_plot[str(16+k*7+2)+"_3D"])<15) & (np.abs(df_plot[str(16+k*7+2)+"_prem_3D"])<15)]

for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+k*7+2)+df])<12)]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*7+1)+df])>0.75)]
    df_plot=df_plot[(np.abs(np.log(df_plot[str(16+k*7+5)+df]))<1.5)]
    df_plot=df_plot[(np.abs(np.log(df_plot[str(16+k*7+3)+df]))<1.5)]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*7+4)+df])<1.5)]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*7+6)+df])<1.5)]

df_plot=df_plot[(df_plot["0"].isin(ev_list))]


df_plot = df_plot[(df_plot["7"]>=90) & (df_plot["7"]<=105)]
#print(df_plot)#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
#pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
#results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],2/5.0,3/5.0)  for i,rows in df_plot.iterrows()])
#pool.close()
#pool.join()
results=[]
for i,rows in df_plot.iterrows(): 
    results.append(points(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],2/5.0,3/5.0)) 
#print(results)
start_lat=[]
start_lon=[]
mid_lat=[]
mid_lon=[]
end_lat=[]
end_lon=[]
for r in results:
    if r is not None:
        
        start_lon.append(r[0][0])
        start_lat.append(r[0][1])
        mid_lat.append(r[1][1])
        mid_lon.append(r[1][0])
        end_lon.append(r[2][0])
        end_lat.append(r[2][1])
        
    else:
        start_lat.append([np.nan,np.nan])
        start_lon.append([np.nan,np.nan])
        mid_lat.append(np.nan)
        mid_lon.append(np.nan)
        end_lat.append([np.nan,np.nan])
        end_lon.append([np.nan,np.nan])
        
df_plot['start_lat'] = start_lat
df_plot['start_lon'] = start_lon
df_plot['mid_lon']   = mid_lon
df_plot['mid_lat']   = mid_lat
df_plot['end_lat']   = end_lat
df_plot['end_lon']   = end_lon

for i,df in enumerate(file_names):

        

    ### Plotting Global Maps for Amplitude Ratio
    fig=pygmt.Figure()
    pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
    fig.basemap(region="d",projection="N-150/12c",frame=True)
    
    
    fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
    fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
    fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(16+k*7+5)+df]),pen="0.000001p",cmap=True,transparency=40)
    fig.coast(shorelines=True, frame=True)
    if i==0:
        fig.colorbar()
    fig.savefig(plot_dir+"/global_amp_"+file_names[i]+"_"+ph+"_high.png",dpi=600)
    plt.close()


plt.figure(1,figsize=(7.08,7.08))


keys = list(zip(df_plot["5"], df_plot["6"]))  # Create a list of tuples where each tuple is (5th column value, 6th column value)  # Get the keys from column "5" where column "6" is True
print(keys)
values = [moho[math.floor(key[1])][math.floor(key[0])] for key in keys]  # Get the corresponding values from the moho dictionary
plt.hist(values, bins=11)#, edgecolor=colors[i], linewidth=lss[i], linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])

#plt.xlim(-1.5,1.5)
plt.legend(fontsize=12)
plt.xlabel("Moho Depth for Stations")
plt.ylabel("Counts")
plt.tight_layout()
plt.savefig(plot_dir+"/Histogram_stations_"+ph+"_high.png",dpi=600)
plt.close()

plt.figure(1,figsize=(7.08,7.08))

    
keys = list(zip(df_plot["3"], df_plot["2"]))  # Create a list of tuples where each tuple is (5th column value, 6th column value)  # Get the keys from column "5" where column "6" is True
values = [moho[math.floor(key[1])][math.floor(key[0])] for key in keys]  # Get the corresponding values from the moho dictionary
plt.hist(values, bins=11,facecolor=None,edgecolor="red",label="Moho Depths")
plt.hist(df_plot["4"],bins=11,facecolor=None,edgecolor="blue",label="Source Depth")

#plt.xlim(-1.5,1.5)
plt.legend(fontsize=12)
plt.xlabel("Moho depth  of Sources")
plt.ylabel("Counts")
plt.tight_layout()
plt.savefig(plot_dir+"/Histogram_events_"+ph+"_high.png",dpi=600)
plt.close()

