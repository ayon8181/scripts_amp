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
#import geopy.distance
from geographiclib.geodesic import Geodesic
from matplotlib.lines import Line2D
from mpi4py import MPI

file_names      = ["_real","_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['darkred','cyan','yellowgreen','chocolate','black']    
lss   =[2,2.5,2,1.5,1] 
style =['None','solid',"solid","solid","solid"]
labels=["real_data",'S40RTS_1D_Crust','S40RTS_3D_crust',"GLAD_M25","PREM_Crust2.0"] #"PREM_QRFSI12",
fccs  =['hotpink','None','None','None',"None"] 
markers=[".","+","o","*","x"]
alp   =[0.5,1,1,1,1]
df=[]

plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/plots_filter"

ev_list=[]
with open("./same.txt","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])

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

phase_list=['S','SS','SSS','ScS','Sdiff']
for k,ph in enumerate(phase_list):
    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all.txt",skipinitialspace=True,delimiter=",")
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

    if k==0:    
        df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
    elif k==1:
        df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
    elif k==2:
        df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<150))]
    elif k==3:
        df_plot = df_plot[((df_plot["7"]>5) & (df_plot["7"]<25)) | ((df_plot["7"]>50) & (df_plot["7"]<65))]
    elif k==4:
        df_plot = df_plot[(df_plot["7"]>100) & (df_plot["7"]<150)]

    ls=["7"]
    df_plot_2=df_plot.copy()
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        ls.append(str(16+k*7+5)+df)
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
    plt.ylabel("$X_{2}$ for "+ph+" phase")
    plt.tight_layout()
    plt.savefig(plot_dir+"/ep_r_"+ph+".png",dpi=600)
    plt.close()

    df_plot_2=df_plot.copy()
    ls=["7"]
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        ls.append(str(16+k*7+3)+df)
    df_dist = df_plot_2[ls].copy()
    for i,cl in enumerate(ls[1:]):
        df_dist[cl] = np.log(df_plot_2[cl])
    
    df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
    df_dist_2["dist"] = np.arange(32.5,150,2.5)
    for i,cl in enumerate(ls[1:]):
        plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5)
    plt.ylim(-0.5,0.5)
    plt.xlabel("Epicentral Distance")
    plt.ylabel("Energy Ratio for "+ph+" phase")
    plt.tight_layout()
    plt.savefig(plot_dir+"/ep_a_"+ph+".png",dpi=600)
    plt.close()
    
    df_plot_2=df_plot.copy()
    ls=["7"]
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        ls.append(str(16+k*7+4)+df)
    df_dist = df_plot_2[ls].copy()
    for cl in ls[1:]:
        df_dist[cl] = df_plot_2[cl]
    df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
    df_dist_2["dist"] = np.arange(32.5,150,2.5)
    for i,cl in enumerate(ls[1:]):
        plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5)
    plt.ylim(-0.5,0.5)
    plt.xlabel("Epicentral Distance")
    plt.ylabel("$X_{env}$ for "+ph+" phase")
    plt.tight_layout()
    plt.savefig(plot_dir+"/ep_e_"+ph+".png",dpi=600)
    plt.close()

    df_plot_2=df_plot.copy()
    ls=["7"]
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        ls.append(str(16+k*7+6)+df)
    df_dist = df_plot_2[ls].copy()
    for cl in ls[1:]:
        df_dist[cl] = df_plot_2[cl]
    df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
    df_dist_2["dist"] = np.arange(32.5,150,2.5)
    for i,cl in enumerate(ls[1:]):
        plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5)
    plt.ylim(-0.5,0.5)
    plt.xlabel("Epicentral Distance")
    plt.ylabel("$X_{1}$ for "+ph+" phase")
    plt.tight_layout()
    plt.savefig(plot_dir+"/ep_am_"+ph+".png",dpi=600)
    plt.close()