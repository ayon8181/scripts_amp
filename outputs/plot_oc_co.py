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
import seaborn as sns
import scipy.stats as stats
font = { "family":"serif",
          "color": "darkred",
          "weight":"normal",
          }

SMALL_SIZE = 18
MEDIUM_SIZE = 24
BIGGER_SIZE = 50

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)

ev_list=[]
with open("./../event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])



file_names      = ["_real","_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['hotpink','cyan','brown','black','yellowgreen']    
lss   =[2,4,3,2,1] 
style =['None','solid',"solid","solid","dashdot"]
labels=["real_data",'S40RTS_1D_Crust','S40RTS_3D_crust',"GLAD_M25","PREM_Crust2.0"] #"PREM_QRFSI12",
fccs  =['hotpink','None','None','None',"None"] 
df=[]
plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/plots_oc"
model=TauPyModel(model="prem")
moho        = {}
with open("/scratch1/09038/ayon8181/scripts_amp/seis/depthtomoho.xyz",'r') as txt:
    data = csv.reader(txt, skipinitialspace=True, delimiter=" ")
    for row in data:
        lats   = math.floor(float(row[1]))
        lons   = math.floor(float(row[0]))
        if lats not in moho.keys():
            moho[lats] = {}
        moho[lats][lons] = float(row[2])


def oc_co(evla, evlo, evdp, stla, stlo, phase_list, moho=moho, model=model):
    #print("run")
    out = np.zeros(len(phase_list))
    out[:] = np.nan
    for j,p in enumerate(phase_list):
        k = model.get_ray_paths_geo(source_depth_in_km=evdp, source_latitude_in_deg=evla, source_longitude_in_deg=evlo, receiver_latitude_in_deg=stla, receiver_longitude_in_deg=stlo,phase_list=[p],resample=True) 
        d3=0
        d1=0
        if len(k) != 0:
            path=k[0].path
            for i,pts in enumerate(path[:-2]):
                if pts[3]==0.0:
                   lat=math.floor(pts[4])
                   lon=math.floor(pts[5])
                   if moho[lat][lon] >= -24.4:
                      d3+=1

            out[j]=d3    
    return out 
phase_list=["SS","SSS"]
for k,ph in enumerate(phase_list):
    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all.txt",skipinitialspace=True,delimiter=",")
    

    
    
    
    
    if k==0:
           df_plot=df[df[str(16+1*6+5)+"_real"].notna()]
           df_plot=df_plot[(df_plot["0"].isin(ev_list))]
           for i,df in enumerate(file_names):
               df_plot=df_plot[(np.abs(df_plot[str(16+1*6+2)+df])<15)]
               df_plot=df_plot[(np.abs(df_plot[str(16+1*6+1)+df])>0.7)]
           df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(oc_co, [(rows[3],rows[2],rows[4],rows[6],rows[5],[ph])  for i,rows in df_plot.iterrows()])
           pool.close()
           pool.join()

           df_plot[ph] = results
           for k,df in enumerate(file_names):
               fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
               ax.hist(np.log(df_plot[df_plot[ph]==0][str(16+1*6+5)+df]),bins=23,range=(-1.5,1.5),edgecolor="orange",facecolor="orange",label="SS Bounce Point at Continental Crust")
               ax.hist(np.log(df_plot[df_plot[ph]==1][str(16+1*6+5)+df]),bins=23,range=(-1.5,1.5),edgecolor="skyblue",facecolor="None",lw=2,label="SS Bounce Point at Oceanic Crust")
               ax.legend()
               ax.set_xlim(-1.5,1.5)
               ax.set_xlabel("Log of Amplitude Ratio")
               ax.set_ylabel("Count")
               plt.savefig(plot_dir+"/Histogram_occo_"+ph+"_"+df+".png",dpi=600)
               plt.close()
    elif k==1:
           df_plot=df[df[str(16+2*6+5)+"_real"].notna()]
           df_plot=df_plot[(df_plot["0"].isin(ev_list))]
           for i,df in enumerate(file_names):
               df_plot=df_plot[(np.abs(df_plot[str(16+2*6+2)+df])<15)]
               df_plot=df_plot[(np.abs(df_plot[str(16+2*6+1)+df])>0.7)]
           df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<165))]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(oc_co, [(rows[3],rows[2],rows[4],rows[6],rows[5],[ph])  for i,rows in df_plot.iterrows()])
           pool.close()
           pool.join()
           df_plot[ph]=results
           for k,df in enumerate(file_names):
               fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
               ax.hist(np.log(df_plot[df_plot[ph]==0][str(16+2*6+5)+df]),bins=23,range=(-1.5,1.5),edgecolor="orange",facecolor="orange",label="Both SSS Bounce Points at Continental Crust")
               ax.hist(np.log(df_plot[df_plot[ph]==1][str(16+2*6+5)+df]),bins=23,range=(-1.5,1.5),edgecolor="black",facecolor="None",lw=4,label="SSS Bounce Point at Continetal and Oceanic Crust")
               ax.hist(np.log(df_plot[df_plot[ph]==2][str(16+2*6+5)+df]),bins=23,range=(-1.5,1.5),edgecolor="skyblue",facecolor="None",lw=2,label="Both SSS Bounce Points at Oceanic Crust")
               ax.set_xlim(-1.5,1.5)
               ax.legend()
               ax.set_xlabel("Log of Amplitude Ratio")
               ax.set_ylabel("Count")
               plt.savefig(plot_dir+"/Histogram_occo_"+ph+"_"+df+".png",dpi=600)
               plt.close()
    

df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_plot=df[df[str(8)+"_real"].notna()]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]
df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
#df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+0*6+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+0*6+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(16+1*6+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+1*6+1)+df])>0.7)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(oc_co, [(rows[3],rows[2],rows[4],rows[6],rows[5],["SS"])  for i,rows in df_plot.iterrows()])
pool.close()
pool.join()
df_plot["SS"]=results

for k,df in enumerate(file_names):
    sns.set_theme(style="darkgrid")
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
    ax.hist(np.log(df_plot[df_plot["SS"]==0]["8"+df]),bins=23,range=(-1.5,1.5),edgecolor="orange",facecolor="orange",label="SS Bounce Point at Continental Crust")
    ax.hist(np.log(df_plot[df_plot["SS"]==1]["8"+df]),bins=23,range=(-1.5,1.5),edgecolor="skyblue",facecolor="None",lw=2,label="SS Bounce Point at Oceanic Crust")
    ax.legend()
    ax.set_xlim(-1.5,1.5)
    ax.set_xlabel("Log of Amplitude Ratio")
    ax.set_ylabel("Count")
    plt.savefig(plot_dir+"/Histogram_occo_SS_S_"+df+".png",dpi=600)
    plt.close()

df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS_SS_all.txt",skipinitialspace=True,delimiter=",")
df_plot=df[df[str(10)+"_real"].notna()]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]
#df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<165))]
for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+2*6+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+2*6+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(16+1*6+2)+df])<15)]
    df_plot=df_plot[(np.abs(df_plot[str(16+1*6+1)+df])>0.7)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(oc_co, [(rows[3],rows[2],rows[4],rows[6],rows[5],["SS"])  for i,rows in df_plot.iterrows()])
pool.close()
pool.join()
df_plot["SS"]=results

pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(oc_co, [(rows[3],rows[2],rows[4],rows[6],rows[5],["SSS"])  for i,rows in df_plot.iterrows()])
pool.close()
pool.join()
df_plot["SSS"]=results

for k,df in enumerate(file_names):
    sns.set_theme(style="darkgrid")
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
    ax.hist(np.log(df_plot[df_plot["SSS"]==0]["10"+df]),bins=23,range=(-1.5,1.5),edgecolor="orange",facecolor="orange",label="Both SSS Bounce Points at Continental Crust")
    ax.hist(np.log(df_plot[df_plot["SSS"]==1]["10"+df]),bins=23,range=(-1.5,1.5),edgecolor="seagreen",facecolor="None",lw=4,label="SSS Bounce Point at Continetal and Oceanic Crust")
    ax.hist(np.log(df_plot[df_plot["SSS"]==2]["10"+df]),bins=23,range=(-1.5,1.5),edgecolor="skyblue",facecolor="None",lw=2,label="Both SSS Bounce Points at Oceanic Crust")
    ax.legend()
    ax.set_xlim(-1.5,1.5)
    ax.set_xlabel("Log of Amplitude Ratio")
    ax.set_ylabel("Count")
    plt.savefig(plot_dir+"/Histogram_occo_SSS_SS_1"+df+".png",dpi=600)
    plt.close()

for k,df in enumerate(file_names):
    sns.set_theme(style="darkgrid")
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
    ax.hist(np.log(df_plot[df_plot["SS"]==0]["10"+df]),bins=23,range=(-1.5,1.5),edgecolor="orange",facecolor="orange",label="SS Bounce Point at Continental Crust")
    ax.hist(np.log(df_plot[df_plot["SS"]==1]["10"+df]),bins=23,range=(-1.5,1.5),edgecolor="skyblue",facecolor="None",lw=2,label="SS Bounce Point at Oceanic Crust")
    ax.legend()
    ax.set_xlim(-1.5,1.5)
    ax.set_xlabel("Log of Amplitude Ratio")
    ax.set_ylabel("Count")
    plt.savefig(plot_dir+"/Histogram_occo_SSS_SS_2"+df+".png",dpi=600)
    plt.close()



