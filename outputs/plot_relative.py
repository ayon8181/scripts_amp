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

font = { "family":"serif",
          "color": "darkred",
          "weight":"normal",
          }

SMALL_SIZE = 16
MEDIUM_SIZE = 24
BIGGER_SIZE = 50

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

ev_list=[]
with open("./../event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])


    
         

file_names      = ["_real","_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['darkred','cyan','yellowgreen','chocolate','black']    
lss   =[2,2.5,2,1.5,1] 
style =['None','solid',"solid","solid","solid"]
labels=["real_data",'S40RTS_1D_Crust','S40RTS_3D_crust',"GLAD_M25","PREM_Crust2.0"] #"PREM_QRFSI12",
fccs  =['hotpink','None','None','None',"None"] 
markers=[".","+","o","*","x"]
alp   =[0.5,1,1,1,1]
df=[]
plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/plots_relative/all"
model=TauPyModel(model="prem")
#for i,f in enumerate(file_names):
#    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt",header=None,delimiter=" ",skipinitialspace=True)
 #   df.append(temp)
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

handles=[]
for i,df in enumerate(file_names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker=markers[i], label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends")
plt.close()


df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_plot=df[df[str(8)+"_real"].notna()]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]

df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
#df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+0*7+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+0*7+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(16+1*7+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+1*7+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(np.log((df_plot[str(8)+df])))<1.5)]
#df_plot=df_plot[(np.abs(df_plot[str(18)+"_real"])<15) & (np.abs(df_plot[str(18)+"_1D"])<15) &(np.abs(df_plot[str(18)+"_3D"])<15) & (np.abs(df_plot[str(18)+"_prem_3D"])<15)]
#df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<15) & (np.abs(df_plot[str(24)+"_1D"])<15) &(np.abs(df_plot[str(24)+"_3D"])<15) & (np.abs(df_plot[str(24)+"_prem_3D"])<15)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],2/5.0,3/5.0)  for i,rows in df_plot.iterrows()])
pool.close()
pool.join()
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

ls=["7"]
df_plot_2=df_plot.copy()
plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    ls.append(str(8)+df)
print(ls)
df_dist = df_plot_2[ls].copy()
for cl in ls[1:]:
    df_dist[cl] = np.log(df_plot_2[cl])
df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
df_dist_2["dist"] = np.arange(32.5,150,2.5)
print(df_dist)
for i,cl in enumerate(ls[1:]):
    plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5)
plt.ylim(-0.5,0.5)
plt.xlabel("Epicentral Distance")
plt.ylabel("SS/S ratio")
plt.tight_layout()
plt.savefig(plot_dir+"/ep_SS_S.png",dpi=600)
plt.close()

for i,df in enumerate(file_names):

        

        ### Plotting Global Maps for Amplitude Ratio
    fig=pygmt.Figure()
    pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
    fig.basemap(region="d",projection="N-150/12c",frame=True)
    
    
    fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
    fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
    fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(8)+df]),pen="0.000001p",cmap=True,transparency=40)
    fig.coast(shorelines=True, frame=True)
    if i==0:
        fig.colorbar()
    fig.savefig(plot_dir+"/global_amp_"+file_names[i]+"_SS_S.png",dpi=600)
    plt.close()

plt.figure(1,figsize=(7.08,7.08))
for i,df in enumerate(file_names):
    
    
    plt.hist(np.log(df_plot[str(8)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i], linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])

plt.xlim(-1.5,1.5)
plt.legend(fontsize=12)
plt.xlabel("SS/S Amplitude Ratio")
plt.ylabel("Counts")
plt.tight_layout()
plt.savefig(plot_dir+"/Histogram_SS_S.png",dpi=600)
plt.close()

plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    
    y, edges = np.histogram(np.log(df_plot[str(8)+df]), bins=23, range=(-1.5,1.5))
    centers = 0.5 * (edges[1:] + edges[:-1])
    plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i])

plt.text(1.3,0.8*max(y/np.sum(y)),"SS/S")
plt.xlim(-1.5,1.5)
#plt.legend(fontsize=12)
plt.xlabel("SS/S amplitude ratio")
plt.ylabel("Normalized Counts")
plt.tight_layout()
plt.grid()
plt.savefig(plot_dir+"/line_SS_S.png",dpi=600)
plt.close()

for i,df in enumerate(file_names[1:]):
       
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
    
    ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
    ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
    ax.scatter(np.log(df_plot[str(8)+"_real"]),np.log(df_plot[str(8)+df]),marker=".",alpha=0.6)
    ax.set_xlabel("SS/S for Real Data")
    ax.set_ylabel("SS/S for "+labels[i+1])
    
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    plt.tight_layout()
    plt.savefig(plot_dir+"/Scatter_SS_S"+"_"+labels[i+1]+".png",dpi=600)
    plt.close()

df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS_SS_all.txt",skipinitialspace=True,delimiter=",")
df_plot=df[df[str(10)+"_real"].notna()]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]

df_plot = df_plot[(((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<140)))]
#df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+1*7+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+1*7+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(16+2*7+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+2*7+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(np.log((df_plot[str(10)+df])))<1.5)]
#df_plot=df_plot[(np.abs(df_plot[str(18)+"_real"])<15) & (np.abs(df_plot[str(18)+"_1D"])<15) &(np.abs(df_plot[str(18)+"_3D"])<15) & (np.abs(df_plot[str(18)+"_prem_3D"])<15)]
#df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<15) & (np.abs(df_plot[str(24)+"_1D"])<15) &(np.abs(df_plot[str(24)+"_3D"])<15) & (np.abs(df_plot[str(24)+"_prem_3D"])<15)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],2/5.0,3/5.0)  for i,rows in df_plot.iterrows()])
pool.close()
pool.join()
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

ls=["7"]
df_plot_2=df_plot.copy()
plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    ls.append(str(10)+df)
print(ls)
df_dist = df_plot_2[ls].copy()
for cl in ls[1:]:
    df_dist[cl] = np.log(df_plot_2[cl])
df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
df_dist_2["dist"] = np.arange(32.5,150,2.5)
print(df_dist)
for i,cl in enumerate(ls[1:]):
    plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5)
plt.ylim(-0.5,0.5)
plt.xlabel("Epicentral Distance")
plt.ylabel("SSS/SS ratio")
plt.tight_layout()
plt.savefig(plot_dir+"/ep_SSS_SS.png",dpi=600)
plt.close()

for i,df in enumerate(file_names):

        

        ### Plotting Global Maps for Amplitude Ratio
    fig=pygmt.Figure()
    pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
    fig.basemap(region="d",projection="N-150/12c",frame=True)
    
    
    fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
    fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
    fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(10)+df]),pen="0.000001p",cmap=True,transparency=40)
    fig.coast(shorelines=True, frame=True)
    if i==0:
        fig.colorbar()
    fig.savefig(plot_dir+"/global_amp_"+file_names[i]+"_SSS_SS.png",dpi=600)
    plt.close()

plt.figure(1,figsize=(7.08,7.08))
for i,df in enumerate(file_names):
    
    
    plt.hist(np.log(df_plot[str(10)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i], linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])

plt.xlim(-1.5,1.5)
plt.legend(fontsize=12)
plt.xlabel("SSS/SS Amplitude Ratio")
plt.ylabel("Counts")
plt.tight_layout()
plt.savefig(plot_dir+"/Histogram_SSS_SS.png",dpi=600)
plt.close()

plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    
    y, edges = np.histogram(np.log(df_plot[str(10)+df]), bins=23, range=(-1.5,1.5))
    centers = 0.5 * (edges[1:] + edges[:-1])
    plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i])

plt.text(1.3,0.8*max(y/np.sum(y)),"SSS/SS")
plt.xlim(-1.5,1.5)
#plt.legend(fontsize=12)
plt.xlabel("SSS/SS amplitude ratio")
plt.ylabel("Normalized Counts")
plt.tight_layout()
plt.grid()
plt.savefig(plot_dir+"/line_SSS_SS.png",dpi=600)
plt.close()

for i,df in enumerate(file_names[1:]):
       
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
    
    ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
    ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
    ax.scatter(np.log(df_plot[str(10)+"_real"]),np.log(df_plot[str(10)+df]),marker=".",alpha=0.6)
    ax.set_xlabel("SSS/SS for Real Data")
    ax.set_ylabel("SSS/SS for "+labels[i+1])
    
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    plt.tight_layout()
    plt.savefig(plot_dir+"/Scatter_SSS_SS"+"_"+labels[i+1]+".png",dpi=600)
    plt.close()
       