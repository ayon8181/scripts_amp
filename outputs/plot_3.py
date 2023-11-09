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

file_names      = ["obsd_1D_crust", "obsd_3D_crust","obsd_glad","real_data"]
df=[]
plot_dir="/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/plots_deep"

for i,f in enumerate(file_names):
    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt2",header=None,delimiter=" ",skipinitialspace=True)
    temp=temp[temp[4]>=100]
    df.append(temp)

ocean_depth = -24
max_cc        = [0.85,0.75,0.75,0.7]
max_tshift    = 15
moho        = {}
with open("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/depthtomoho.xyz",'r') as txt:
    data = csv.reader(txt, skipinitialspace=True, delimiter=" ")
    for row in data:
        lats   = math.floor(float(row[1]))
        lons   = math.floor(float(row[0]))
        if lats not in moho.keys():
            moho[lats] = {}
        moho[lats][lons] = float(row[2])   
        

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

model=TauPyModel(model="prem")
phase_list=["SS","SSS"]
def oc_co(evla, evlo, evdp, stla, stlo, moho=moho, model=model, phase_list=phase_list):
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


for d in df:
    pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
    results=pool.starmap(oc_co, [(rows[3],rows[2],rows[4],rows[6],rows[5])  for i,rows in d.iterrows()])
    pool.close()
    pool.join()
    results=np.reshape(results,(len(results),2))
    for i,p in enumerate(phase_list):
        d[p] = results[:,i]
    print(d)

colors=['green','red','black','blueviolet']    
lss   =[3,3,3,3] 
labels=['S40RTS_1D_Crust', 'S40RTS_3D_crust','GLAD_M25',"real_data"] 
fccs  =['None','None','None','none'] 
phase_list=["S","SS","SSS","ScS","Sdiff"]



plt.figure(1, figsize=[10,10])
for i,dfs in enumerate(df):
    df_plot = dfs[(dfs[18]>max_cc[i]) & (dfs[23]>max_cc[i]) & (np.abs(dfs[19]<15)) & (np.abs(dfs[24]<15))]
    df_plot = df_plot[df_plot[8] == df_plot[8]]
    #print(df_plot[14])
    plt.hist(np.log(df_plot[8]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i],alpha=0.5)
    plt.axvline(np.mean(np.log(df_plot[8])),color=colors[i],linestyle="dashed",linewidth=2,alpha=0.5)
plt.legend()
plt.xlabel("Log of SS/S amplitude ratio")
plt.ylabel("Counts")
plt.savefig(plot_dir+"/Histogram_SS_S.png")
plt.close()

plt.figure(1,figsize=[10,10])
for i,dfs in enumerate(df):
    df_plot = dfs[(dfs[23]>max_cc[i]) & (dfs[28]>max_cc[i]) & (np.abs(dfs[24]<15)) & (np.abs(dfs[29]<15))]
    df_plot = df_plot[df_plot[10] == df_plot[10]]
    plt.hist(np.log(df_plot[10]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i],alpha=0.5)
    plt.axvline(np.mean(np.log(df_plot[10])),color=colors[i],linestyle="dashed",linewidth=2)
plt.legend()
plt.xlabel("Log of SSS/SS amplitude ratio")
plt.ylabel("Counts")
plt.savefig(plot_dir+"/Histogram_SSS_SS.png")
plt.close()

for k,ph in enumerate(phase_list[1:3]):
    
    for i,dfs in enumerate(df):
        

        #if k==0:
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[ph] == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        if k==0:
           df_plot = dfs[(dfs[18]>max_cc[i]) & (dfs[23]>max_cc[i]) & (np.abs(dfs[19]<15)) & (np.abs(dfs[24]<15))]
           df_plot = df_plot[df_plot[8] == df_plot[8]]
           plt.figure(1, figsize=[10,10])
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<70)][8]), bins=23, range=(-1.5,1.5), edgecolor='blueviolet', linewidth=3, label="All Paths",facecolor='None')
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<70) & (df_plot['SS']  == 0)][8]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=2, label="SS bounce point at Continental crust",facecolor='None', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<70) & (df_plot['SS']  == 1)][8]), bins=23, range=(-1.5,1.5), edgecolor='green', linewidth=1, label="SS bounce points at Oceanic crust",facecolor='None', alpha = 0.5)
           plt.legend()
           plt.xlabel("Log of SS/S amplitude Ratios for"+labels[i])
           plt.ylabel("Counts")
           plt.savefig(plot_dir+"/Histogram_SS_S_"+labels[i]+".png")
           plt.close()
        elif k==1:
           df_plot = dfs[(dfs[23]>max_cc[i]) & (dfs[28]>max_cc[i]) & (np.abs(dfs[24]<15)) & (np.abs(dfs[29]<15))]
           df_plot = df_plot[df_plot[10] == df_plot[10]]
           plt.figure(1, figsize=[10,10])
           plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140))][10]), bins=23, range=(-1.5,1.5), edgecolor='blueviolet', linewidth=4,label="All paths",facecolor='None')
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140))) & (df_plot['SSS']  == 0)][10]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=3, label="Both SSS bounce points at Continental Crusts",facecolor='None', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140))) & (df_plot['SSS']  == 1)][10]), bins=23, range=(-1.5,1.5), edgecolor='green', linewidth=2, label="One SSS bounce points at Continental and other at Oceanic Crust ",facecolor='None', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140))) & (df_plot['SSS']  == 2)][10]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1,label="Both SSS bounce points at Oceanic Crusts",facecolor='None', alpha = 0.5)
           plt.legend()
           plt.xlabel("Log of SSS/SS amplitude Ratios for"+labels[i])
           plt.ylabel("Counts")
           plt.savefig(plot_dir+"/Histogram_1_SSS_SS_"+labels[i]+".png")
           plt.close()
           plt.figure(1, figsize=[10,10])
           plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140))][10]), bins=23, range=(-1.5,1.5), edgecolor='blueviolet', linewidth=3, label="All Paths",facecolor='blueviolet')
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140))) & (df_plot['SS']  == 0)][10]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=2, label="SS bounce point at Continental crust",facecolor='None', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140))) & (df_plot['SS']  == 1)][10]), bins=23, range=(-1.5,1.5), edgecolor='green', linewidth=1, label="SS bounce points at Oceanic crust",facecolor='None', alpha = 0.5)
           plt.legend()
           plt.xlabel("Log of SSS/SS amplitude Ratios for"+labels[i])
           plt.ylabel("Counts")
           plt.savefig(plot_dir+"/Histogram_2_SSS_SS_"+labels[i]+".png")
           plt.close()

#df=df[:-2]


for k,ph in enumerate(phase_list):
    plt.figure(1,figsize=[10,10])
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[16+k*6+1]>max_cc[i]) & (np.abs((dfs[16+k*6+2])<15))]
        df_plot = df_plot[df_plot[16+k*6+3] == df_plot[16+k*6+3]]
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+3])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==1:
             plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+3])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==2:
            plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
            plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+3])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==3:
             plt.hist(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+3])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==4:
             plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5) 
             plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+3])),color=colors[i],linestyle="dashed",linewidth=2)
    plt.legend()
    plt.xlabel("Log of 3d/1d amplitude ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_3d_1d_"+ph+".png")
    plt.close()

for k,ph in enumerate(phase_list):
    plt.figure(1, figsize=[10,10])
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[16+k*6+1]>max_cc[i]) & (np.abs((dfs[16+k*6+2])<15))]
        df_plot = df_plot[df_plot[16+k*6+4] == df_plot[16+k*6+4]]
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+4]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+4])),color=colors[i],linestyle="dashed",linewidth=2)      
        elif k==1:
             plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+4]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+4])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==2:
             plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165)) ][16+k*6+4]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+4])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==3:
             plt.hist(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+4]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+4])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==4:
             
             plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+4]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)    
             plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+4])),color=colors[i],linestyle="dashed",linewidth=2)
    plt.legend()
    plt.xlabel("Log of 3d/1d Envelope ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_envp_3d_1d_"+ph+".png")
    plt.close()

for k,ph in enumerate(phase_list):
    plt.figure(1, figsize=[10,10])
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[16+k*6+1]>max_cc[i]) & (np.abs((dfs[16+k*6+2])<15))]
        df_plot = df_plot[df_plot[16+k*6+2] == df_plot[16+k*6+2]]
        if k==0:
           plt.hist(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+2], bins=23, range=(-15,15), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
           plt.axvline(np.mean(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+2]),color=colors[i],linestyle="dashed",linewidth=2)      
        elif k==1:
             plt.hist(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+2], bins=23, range=(-15,15), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+2]),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==2:
             plt.hist(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165)) ][16+k*6+2], bins=23, range=(-15,15), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+2]),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==3:
             plt.hist(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+2], bins=23, range=(-15,15), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))][16+k*6+2]),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==4:
             
             plt.hist(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+2], bins=23, range=(-15,15), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)    
             plt.axvline(np.mean(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+2]),color=colors[i],linestyle="dashed",linewidth=2)
    plt.legend()
    plt.xlabel("3D-1D CC travel time difference "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_shift_"+ph+".png")
    plt.close()

for k,ph in enumerate(phase_list):
    plt.figure(1, figsize=[10,10])
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[16+k*6+1]>max_cc[i]) & (np.abs((dfs[16+k*6+2])<15))]
        df_plot = df_plot[df_plot[16+k*6+5] == df_plot[16+k*6+5]]
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+5])),color=colors[i],linestyle="dashed",linewidth=2)      
        elif k==1:
             plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+5])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==2:
             plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165)) ][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+5])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==3:
             plt.hist(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
             plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+5])),color=colors[i],linestyle="dashed",linewidth=2)
        elif k==4:
             
             plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)    
             plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150)][16+k*6+4])),color=colors[i],linestyle="dashed",linewidth=2)
    plt.legend()
    plt.xlabel("Log of 3d/1d Amplitude Ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_amp_rat"+ph+".png")
    plt.close()
"""

for k,ph in enumerate(phase_list):
    
    for i,dfs in enumerate(df):
        plt.figure(1, figsize=[10,10])
        df_plot = dfs[(dfs[16+k*6+2]>0.85) & (np.abs(dfs[16+k*6+3]<15))]
        df_plot = df_plot[df_plot[16+k*6+4] == df_plot[16+k*6+4]]
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[13]>ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[13]<=ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='blueviolet', alpha = 0.5)
        elif k==1:
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[13]>ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[13]<=ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='blueviolet', alpha = 0.5)
        elif k==2:
           plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[13]>ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[13]<=ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='blueviolet', alpha = 0.5)
        elif k==3:
           plt.hist(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[13]>ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[13]<=ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='blueviolet', alpha = 0.5)
        elif k==4:
           plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[13]>ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[13]<=ocean_depth)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='blueviolet', alpha = 0.5)
        plt.legend()
        plt.xlabel("Log of 3d/1d Amplitude ratio for "+ph+" phase for "+labels[i])
        plt.ylabel("Counts")
        plt.savefig(plot_dir+"/Histogram_"+labels[i]+"_3d_1d_"+ph+".png")
        plt.close()

"""
for k,ph in enumerate(phase_list):
    
    for i,dfs in enumerate(df):
        plt.figure(1, figsize=[10,10])
        df_plot = dfs[(dfs[16+k*6+1]>max_cc[i]) & (np.abs(dfs[16+k*6+2]<15))]
        df_plot = df_plot[df_plot[16+k*6+3] == df_plot[16+k*6+3]]
        #if k==0:
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[ph] == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        if k==1:
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blueviolet', linewidth=3,label="All Paths",facecolor='None')
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+3])),color="blueviolet",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=2,label="SS bounce point at Continental crust",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 0)][16+k*6+3])),color="red",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='green', linewidth=1,label="SS bounce points at Oceanic crust",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 1)][16+k*6+3])),color="green",linestyle="dashed",linewidth=1)
        elif k==2:
           plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blueviolet', linewidth=4,label="All Paths",facecolor='None')
           plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+3])),color="blueviolet",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=3,label="Both SSS bounce points at Continental Crusts",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 0)][16+k*6+3])),color="red",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='green', linewidth=2, label="One SSS bounce points at Continental and other at Oceanic Crust ",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 1)][16+k*6+3])), color="green",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 2)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1,label="Both SSS bounce points at Oceanic Crusts",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 2)][16+k*6+3])), color="black",linestyle="dashed",linewidth=1)
        #elif k==3:
        #   plt.hist(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        #   plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        #elif k==4:
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        plt.legend()
        plt.xlabel("Log of 3d/1d Amplitude ratio for "+ph+" phase for "+labels[i])
        plt.ylabel("Counts")
        plt.savefig(plot_dir+"/Histogram_"+labels[i]+"_3d_1d_"+ph+".png")
        plt.close()

for k,ph in enumerate(phase_list):
    
    for i,dfs in enumerate(df):
        plt.figure(1, figsize=[10,10])
        df_plot = dfs[(dfs[16+k*6+1]>max_cc[i]) & (np.abs(dfs[16+k*6+2]<15))]
        df_plot = df_plot[df_plot[16+k*6+5] == df_plot[16+k*6+5]]
        #if k==0:
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[ph] == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>30) & (df_plot[7]<70) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        if k==1:
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor='blueviolet', linewidth=3,label="All Paths",facecolor='None')
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140)][16+k*6+5])),color="blueviolet",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 0)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=2,label="SS bounce point at Continental crust",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 0)][16+k*6+5])),color="red",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 1)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor='green', linewidth=1,label="SS bounce points at Oceanic crust",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(df_plot[7]>50) & (df_plot[7]<140) & (df_plot[ph]  == 1)][16+k*6+5])),color="green",linestyle="dashed",linewidth=1)
        elif k==2:
           plt.hist(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor='blueviolet', linewidth=4,label="All Paths",facecolor='None')
           plt.axvline(np.mean(np.log(df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))][16+k*6+5])),color="blueviolet",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 0)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=3,label="Both SSS bounce points at Continental Crusts",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 0)][16+k*6+5])),color="red",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 1)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor='green', linewidth=2, label="One SSS bounce points at Continental and other at Oceanic Crust ",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 1)][16+k*6+5])), color="green",linestyle="dashed",linewidth=1)
           plt.hist(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 2)][16+k*6+5]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1,label="Both SSS bounce points at Oceanic Crusts",facecolor='None', alpha = 0.5)
           plt.axvline(np.mean(np.log(df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))) & (df_plot[ph]  == 2)][16+k*6+5])), color="black",linestyle="dashed",linewidth=1)
        #elif k==3:
        #   plt.hist(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        #   plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        #elif k==4:
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        plt.legend()
        plt.xlabel("Log of Amplitude Ratio for "+ph+" phase for "+labels[i])
        plt.ylabel("Counts")
        plt.savefig(plot_dir+"/Histogram_"+labels[i]+"_amp_"+ph+".png")
        plt.close()


        #elif k==3:
        #   plt.hist(np.log(df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        #   plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65))) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)
        #elif k==4:
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) ][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[ph]  == 0)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[7]>100) & (df_plot[7]<150) & (df_plot[ph]  == 1)][16+k*6+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='blueviolet', alpha = 0.5)


