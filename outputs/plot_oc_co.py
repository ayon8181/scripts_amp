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
import mpl_scatter_density

font = { "family":"serif",
          "color": "darkred",
          "weight":"normal",
          }

SMALL_SIZE = 6
MEDIUM_SIZE = 12
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

ev_list=[]
with open("./../event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])


moho={}
with open("/scratch1/09038/ayon8181/scripts_amp/outputs/depthtomoho.xyz",'r') as txt:
    data = csv.reader(txt, skipinitialspace=True, delimiter=" ")
    for row in data:
        lats   = math.floor(float(row[1]))
        lons   = math.floor(float(row[0]))
        if lats not in moho.keys():
            moho[lats] = {}
        moho[lats][lons] = float(row[2])
    moho[90]=moho[89]
    for keys in moho.keys():
        moho[keys][180]=moho[keys][-180]
#print(moho)


us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"]         

file_names      = ["_real","_1D","_3D","_3D_2","_3D_15","_prem_3D","_glad"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['crimson','forestgreen','peru','magenta','indigo','deepskyblue','darkorange','grey','black']     
lss   =[2,1.5,1,0.5,1,1,1,1,1,1,1,1,1,1] 
style =['None','solid',"solid","solid","solid","solid",'solid']
labels=["real_data",'SPREMc','S8000',"S11000","S5000","PREMc","GLAD_M25"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None"] 
markers=[".","+","o","*","x","^","v"]
alp   =[0.5,1,1,1,1,1,1]
df=[]
plot_dir="/scratch1/09038/ayon8181/synthetic/scripts_amp/outputs/figures/oc_co/test"
model=TauPyModel(model="prem")


labels=["PREM_3DC","S40RTS_1DC","S40RTS_3DC","GLAD_M25"]
file_names=["_prem_3D","_1D","_3D","_glad"]

def points(evla,evlo,evdp,stla,stlo,gcarc):
    az=Geodesic.WGS84.Inverse(evla, evlo, stla, stlo)['azi1']
    
    m_di=gcarc*(1/2.0)
    mid=geopy.distance.distance(nautical=60*m_di).destination((evla,evlo),bearing=az)
    
    recv=0
    if moho[math.floor(stla)][math.floor(stlo)]<-24:
        recv=1
    midp=0
    if moho[math.floor(mid[0])][math.floor(mid[1])]<-24:
        midp=1

    return [recv,midp]

handles=[]
for i,df in enumerate(file_names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker=markers[i], label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends")
plt.close()

def plot(ph,k):
    df=pd.read_csv("/scratch1/09038/ayon8181/synthetic/scripts_amp/outputs/"+ph+"_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))]
    df_plot=df_plot[(df_plot["7"]>30)]
    df_plot=df_plot[(df_plot["4"]>50)]
    #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(file_names):
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+df])<20)]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+1)+df])>0.7)]
        df_plot=df_plot[df_plot[str(18+k*11+5)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+4)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+3)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+3)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+4)+df].notna()]
        df_plot=df_plot[(np.abs(np.log(df_plot[str(18+k*11+5)+df])<1.5))]
        df_plot=df_plot[df_plot[str(18+k*11+6)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+6)+df])<1.5)]
        
    df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    #df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]
    print(len(df_plot["1"]))
    df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    (df_plot["1"].str.split(".").str[0].isin(gsn))]
    print(len(df_plot["1"]))

    if k==0:    
        df_plot = df_plot[((df_plot["7"]>30) & (df_plot["7"]<75)) | ((df_plot["7"]>85) & (df_plot["7"]<100))]
        ranges  = [30,60,100]
       
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"])  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        
        recvs=[]
        midps=[]
        for r in results:
            if r is not None:
                
                recvs.append(r[0])
                midps.append(r[1])
                
            else:
                recvs.append(np.nan)
                midps.append(np.nan)
        df_plot['recvs'] = recvs
        df_plot['midps'] = midps
        
    
    elif k==1:
        df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
        ranges  = [50,80,110,140]

        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"])  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        
        recvs=[]
        midps=[]
        for r in results:
            if r is not None:
                
                recvs.append(r[0])
                midps.append(r[1])
                
            else:
                recvs.append(np.nan)
                midps.append(np.nan)
        df_plot['recvs'] = recvs
        df_plot['midps'] = midps
    elif k==2:
        df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>115) & (df_plot["7"]<150))]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        ranges  = [70,100,125,150]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"])  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        
        recvs=[]
        midps=[]
        for r in results:
            if r is not None:
                
                recvs.append(r[0])
                midps.append(r[1])
                
            else:
                recvs.append(np.nan)
                midps.append(np.nan)
        df_plot['recvs'] = recvs
        df_plot['midps'] = midps
    elif k==3:
        df_plot = df_plot[((df_plot["7"]>5) & (df_plot["7"]<25)) | ((df_plot["7"]>50) & (df_plot["7"]<65)) ]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        ranges=[5,65]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"])  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        
        recvs=[]
        midps=[]
        for r in results:
            if r is not None:
                
                recvs.append(r[0])
                midps.append(r[1])
                
            else:
                recvs.append(np.nan)
                midps.append(np.nan)
        df_plot['recvs'] = recvs
        df_plot['midps'] = midps
    elif k==4:
        df_plot = df_plot[(df_plot["7"]>100) & (df_plot["7"]<150)]
        ranges  = [100,125,150]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"])  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        
        recvs=[]
        midps=[]
        for r in results:
            if r is not None:
                
                recvs.append(r[0])
                midps.append(r[1])
                
            else:
                recvs.append(np.nan)
                midps.append(np.nan)
        df_plot['recvs'] = recvs
        df_plot['midps'] = midps

    fig,ax=plt.subplots(2,2,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
        y, edges = np.histogram(df_plot[str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
        #centers = 0.5 * (edges[1:] + edges[:-1])
        if i==3:
            ax[int(i/2),int(i%2)].stairs(y/np.sum(y),edges,color="cornflowerblue",fill="cornflowerblue",label="All")
            y2, edges = np.histogram(df_plot[df_plot["recvs"]==1][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y2/np.sum(y),edges,color="black",label="Continental Crust",lw=2)
            y3, edges = np.histogram(df_plot[df_plot["recvs"]==0][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y3/np.sum(y),edges,color="red",label="Oceanic Crust",lw=1)
            ax[int(i/2),int(i%2)].legend()
            #ax[int(i/2),int(i%2)].text(1.0,0.8*max(y/np.sum(y)),ph)
        else:
            ax[int(i/2),int(i%2)].stairs(y/np.sum(y),edges,color="cornflowerblue",fill="cornflowerblue")
            y2, edges = np.histogram(df_plot[df_plot["recvs"]==1][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y2/np.sum(y),edges,color="black",lw=2)
            y3, edges = np.histogram(df_plot[df_plot["recvs"]==0][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y3/np.sum(y),edges,color="red",lw=1)
            #if i==1:
            #   ax[int(i/2),int(i%2)].text(1.0,0.8*max(y/np.sum(y)),ph)
        ax[int(i/2),int(i%2)].text(0.8,0.3*max(y/np.sum(y)),labels[i])
        ax[int(i/2),int(i%2)].set_xlim(-1.4,1.4)
        #plt.legend(fontsize=12)
        #ax[int(i/2),int(i%2)].set_xlabel("$\\chi_{amp}$")
        #ax[int(i/2),int(i%2)].set_ylabel("Normalized Counts")
        
        ax[int(i/2),int(i%2)].grid()
    fig.supylabel('Normalized Counts of '+ph,fontsize=12)
    fig.supxlabel("$\\chi_{amp}$",fontsize=12)
    plt.savefig(plot_dir+"/receiver_amp_"+ph+".png",dpi=600)
    plt.close()

    fig,ax=plt.subplots(2,2,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
        y, edges = np.histogram(df_plot[str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
        #centers = 0.5 * (edges[1:] + edges[:-1])
        if i==3:
            ax[int(i/2),int(i%2)].stairs(y/np.sum(y),edges,color="cornflowerblue",fill="cornflowerblue",label="All")
            y2, edges = np.histogram(df_plot[df_plot["midps"]==1][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y2/np.sum(y),edges,color="black",label="Continental Crust",lw=2)
            y3, edges = np.histogram(df_plot[df_plot["midps"]==0][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y3/np.sum(y),edges,color="red",label="Oceanic Crust",lw=1)
            ax[int(i/2),int(i%2)].legend()
        
        else:
            ax[int(i/2),int(i%2)].stairs(y/np.sum(y),edges,color="cornflowerblue",fill="cornflowerblue")
            y2, edges = np.histogram(df_plot[df_plot["midps"]==1][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y2/np.sum(y),edges,color="black",lw=2)
            y3, edges = np.histogram(df_plot[df_plot["midps"]==0][str(18+k*11+6)+df], bins=14, range=(-1.4,1.4))
            ax[int(i/2),int(i%2)].stairs(y3/np.sum(y),edges,color="red",lw=1)
            # if i==1:
            #    ax[int(i/2),int(i%2)].text(1.0,0.8*max(y/np.sum(y)),ph)
        ax[int(i/2),int(i%2)].text(0.8,0.3*max(y/np.sum(y)),labels[i])
        #ax[int(i/2),int(i%2)].set_title(labels[i],fontsize=6)
        ax[int(i/2),int(i%2)].set_xlim(-1.4,1.4)
        #plt.legend(fontsize=12)
        
        
        ax[int(i/2),int(i%2)].grid()
    fig.supylabel('Normalized Counts of '+ph,fontsize=12)
    fig.supxlabel("$\\chi_{amp}$",fontsize=12)
    plt.savefig(plot_dir+"/midpoint_amp_"+ph+".png",dpi=600)
    plt.close()

plot('Sdiff',4)
plot('SS',1)
plot('S',0)
plot('SSS',2)
