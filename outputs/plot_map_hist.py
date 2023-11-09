#!/usr/bin/env python
import numpy as np
import csv
import pygmt
import pandas as pd
import os
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
import math

file_names      = ["obsd_1D_crust", "obsd_3D_crust","obsd_glad"]
df=[]
plot_dir="/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/plots"

for i,f in enumerate(file_names):
    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt4",header=None,delimiter=" ",skipinitialspace=True)
    df.append(temp)

ocean_depth = -24
max_cc        = 0.85
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
        


model=TauPyModel(model="prem")
phase_list=['SS','SSS']
def oc_co(evla, evlo, evdp, stla, stlo, moho=moho, model=model, phase_list=phase_list):
    out = np.zeros(len(phase_list))
    out[:] = np.nan
    for j,p in enumerate(phase_list):
        k = model.get_ray_paths_geo(source_depth_in_km=evdp, source_latitude_in_deg=evla, source_longitude_in_deg=evlo, receiver_latitude_in_deg=stla, receiver_longitude_in_deg=stlo,phase_list=[p],resample=True) 
        d3=0
        d1=0
        if len(k) != 0:
            path=k[0].path
            for i,pts in enumerate(path):
                lat=math.floor(pts[4])
                lon=math.floor(pts[5])
                if moho[lat][lon] >= -pts[3]:
                   d3+=pts[2]
                if pts[3] <= 24.4:
                   d1+=pts[2]

            #if d1 > d3:
            #   out[j] = 1             #### For continents it's 1
            #else:
            #     out[j] = 0
            out[j] = d3/d1
    return out 
               
            
    


def percent_ocean(a,o_depth):
    pts = len(a)
    o=0
    c=0
    for j in a:
        if a > o_depth:
            o+=1
        else:
             c+=1
    return (o/pts)



for d in df:
    d['ocean_percent'] = d.apply(lambda x: oc_co(x[2], x[1], x[3], x[5], x[4]), axis=1)

colors=['blue','red','black']    
lss   =['solid','dotted','dashed'] 
labels=['S40RTS_1D_Crust', 'S40RTS_3D_crust','GLAD_M25'] 
fccs  =['skyblue','None','None'] 
phase_list=["S","SS","SSS","ScS","Sdiff"]

plt.figure(1, figsize=[10,10])
for i,dfs in enumerate(df):
    df_plot = dfs[(dfs[24]>0.85) & (dfs[29]>0.85)]
    df_plot = df_plot[df_plot[15] == df_plot[15]]
    #print(df_plot[14])
    plt.hist(np.log(df_plot[15]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i])
plt.legend()
plt.xlabel("Log of SS/S amplitude ratio")
plt.ylabel("Counts")
plt.savefig(plot_dir+"/Histogram_SS_S.png")
plt.close()

plt.figure(1,figsize=[10,10])
for i,dfs in enumerate(df):
    df_plot = dfs[(dfs[34]>0.85) & (dfs[29]>0.85)]
    df_plot = df_plot[df_plot[17] == df_plot[17]]
    plt.hist(np.log(df_plot[17]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i])
plt.legend()
plt.xlabel("Log of SSS/SS amplitude ratio")
plt.ylabel("Counts")
plt.savefig(plot_dir+"/Histogram_SSS_SS.png")
plt.close()



for k,ph in enumerate(phase_list):
    plt.figure(1,figsize=[10,10])
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[23+k*5+1]>0.85) & (np.abs((dfs[23+k*5+2])<15))]
        df_plot = df_plot[df_plot[23+k*5+3] == df_plot[23+k*5+3]]
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==1:
             plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==2:
            plt.hist(np.log(df_plot[((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==3:
             plt.hist(np.log(df_plot[((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65)) ][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==4:
             plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5) 

    plt.legend()
    plt.xlabel("Log of 3d/1d amplitude ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_3d_1d_"+ph+".png")
    plt.close()

    for k,ph in enumerate(phase_list):
        plt.figure(1,figsize=[10,10])
        for i,dfs in enumerate(df):
            df_plot = dfs[(dfs[23+k*5+1]>0.85) & (np.abs((dfs[23+k*5+2])<15))]
            df_plot = df_plot[df_plot[23+k*5+3] == df_plot[23+k*5+3]]
            if k==0:
               plt.hist(df_plot[(df_plot[6]>30) & (df_plot[6]<70)][ph], bins=23, edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
            elif k==1:
                plt.hist(df_plot[(df_plot[6]>50) & (df_plot[6]<140)][ph], bins=23, edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
            elif k==2:
                plt.hist(df_plot[((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))][ph], bins=23, edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
            elif k==3:
                plt.hist(df_plot[((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65)) ][ph], bins=23, edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
            elif k==4:
                plt.hist(df_plot[(df_plot[6]>100) & (df_plot[6]<150)][ph], bins=23, edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5) 

        plt.legend()
        plt.xlabel("Log of 3d/1d amplitude ratio for "+ph+" phase")
        plt.ylabel("Counts")
        plt.savefig(plot_dir+"/Histogram_3d_1d_dist_"+ph+".png")
        plt.close()

for k,ph in enumerate(phase_list):
    plt.figure(1, figsize=[10,10])
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[23+k*5+1]>0.85) & (np.abs((dfs[23+k*5+2])<15))]
        df_plot = df_plot[df_plot[23+k*5+4] == df_plot[23+k*5+4]]
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70)][23+k*5+4]), bins=23, range=(-1,1), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==1:
             plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140)][23+k*5+4]), bins=23, range=(-1,1), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==2:
            plt.hist(np.log(df_plot[((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165)) ][23+k*5+4]), bins=23, range=(-1,1), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==3:
             plt.hist(np.log(df_plot[((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65)) ][23+k*5+4]), bins=23, range=(-1,1), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        elif k==4:
             plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150)][23+k*5+4]), bins=23, range=(-1,1), edgecolor=colors[i], linewidth=1, linestyle=lss[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)    
    plt.legend()
    plt.xlabel("Log of 3d/1d Envelope ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_envp_3d_1d_"+ph+".png")
    plt.close()




for k,ph in enumerate(phase_list[1:3]):
    
    for i,dfs in enumerate(df):
        plt.figure(1, figsize=[10,10])
        df_plot = dfs[(dfs[23+k*5+2]>0.85) & (np.abs(dfs[23+k*5+3]<15))]
        df_plot = df_plot[df_plot[23+k*5+4] == df_plot[23+k*5+4]]
        #if k==0:
        #   plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        #   plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70) & (df_plot[13]>ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
        #   plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70) & (df_plot[13]<=ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='skyblue', alpha = 0.5)
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140) & (df_plot[13]>ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140) & (df_plot[13]<=ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='skyblue', alpha = 0.5)
        elif k==1:
           plt.hist(np.log(df_plot[((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))) & (df_plot[13]>ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))) & (df_plot[13]<=ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='skyblue', alpha = 0.5)
        #elif k==3:
        #   plt.hist(np.log(df_plot[((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65)) ][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
        ##   plt.hist(np.log(df_plot[(((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65))) & (df_plot[13]>ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
         #  plt.hist(np.log(df_plot[(((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65))) & (df_plot[13]<=ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='skyblue', alpha = 0.5)
        ##elif k==4:
          # plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150) ][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None', alpha = 0.5)
          # plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150) & (df_plot[13]>ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="Midpoint at Ocean",facecolor='orange', alpha = 0.5)
           #plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150) & (df_plot[13]<=ocean_depth)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="Midpoint at Continent",facecolor='skyblue', alpha = 0.5)
        plt.legend()
        plt.xlabel("Log of 3d/1d Amplitude ratio for "+ph+" phase for "+labels[i])
        plt.ylabel("Counts")
        plt.savefig(plot_dir+"/Histogram_"+labels[i]+"_3d_1d_"+ph+".png")
        plt.close()


for k,ph in enumerate(phase_list):
    
    for i,dfs in enumerate(df):
        plt.figure(1, figsize=[10,10])
        df_plot = dfs[(dfs[23+k*5+2]>0.85) & (np.abs(dfs[23+k*5+3]<15))]
        df_plot = df_plot[df_plot[23+k*5+4] == df_plot[23+k*5+4]]
        if k==0:
           plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70) & (df_plot['oceanic_depth'][k] == 0)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[6]>30) & (df_plot[6]<70) & (df_plot['oceanic_depth'][k] == 1)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='skyblue', alpha = 0.5)
        elif k==1:
           plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140) & (df_plot['oceanic_depth'][k] == 0)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[6]>50) & (df_plot[6]<140) & (df_plot['oceanic_depth'][k] == 1)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='skyblue', alpha = 0.5)
        elif k==2:
           plt.hist(np.log(df_plot[((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))) & (df_plot['oceanic_depth'][k] == 0)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))) & (df_plot['oceanic_depth'][k] == 1)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='skyblue', alpha = 0.5)
        elif k==3:
           plt.hist(np.log(df_plot[((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65)) ][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None')
           plt.hist(np.log(df_plot[(((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65))) & (df_plot['oceanic_depth'][k] == 0)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65))) & (df_plot['oceanic_depth'][k] == 1)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='skyblue', alpha = 0.5)
        elif k==4:
           plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150) ][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='black', linewidth=1, linestyle=lss[i],label=labels[i],facecolor='None', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150) & (df_plot['oceanic_depth'][k] == 0)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='red', linewidth=1, linestyle='dashed',label="More Oceanic Crust than PREM",facecolor='orange', alpha = 0.5)
           plt.hist(np.log(df_plot[(df_plot[6]>100) & (df_plot[6]<150) & (df_plot['oceanic_depth'][k] == 1)][23+k*5+3]), bins=23, range=(-1.5,1.5), edgecolor='blue', linewidth=1, linestyle='dashed',label="More Continental Crust than PREM",facecolor='skyblue', alpha = 0.5)
        plt.legend()
        plt.xlabel("Log of 3d/1d Amplitude ratio for "+ph+" phase for "+labels[i])
        plt.ylabel("Counts")
        plt.savefig(plot_dir+"/Histogram_paths_"+labels[i]+"_3d_1d_"+ph+".png")
        plt.close()



