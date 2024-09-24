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
#from mpi4py import MPI

us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"]
file_names      = ["_real","_prem_3D","_3D","_3D_2","_3D_15","_1D","_glad","_1D_ref","_prem_q16","_3D_atten","_1D"]
colors=['peru','forestgreen','black','magenta','mediumpurple','deepskyblue','darkorange','grey','crimson']     
lss   =[2,1.5,1.75,1.25,0.5,1,1,1,1,1,1,1,1,1] 
style =['solid','solid',"solid","solid","solid","solid","solid","solid","solid"]
labels=["Real Data","PREMc","S8000","S11000","S5000","SPREMc","GLAD_M25","1D_REF","PREM_QL6.txt","PREM_QRFSI12","S40RTS_PREM_Crust"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None","None","None","None","None"] 
markers=["v",".","+","o","*","x","^"]
alp   =[1,1,1,1,1,1,1]
df=[]

plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/figures"



ev_list=[]
with open("./all.txt","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])






labels=["PREM","PREM_3DC","S40RTS_1DC","S40RTS_3DC","GLAD_M25"]
file_names=["_synt","_prem_3D","_1D","_3D","_glad"]


# labels=["GLAD_M25","S40RTS_3DC","S40RTS_3DC_QRFSI12","PREM","1D_REF","PREM_Q16"]
# file_names=["_glad","_3D","_S80_3D","_synt","_1D_ref","_prem_16"]
all_f=["_synt","_prem_3D","_1D","_3D","_glad"]


#all_f=["_glad","_3D","_S80_3D","_synt","_1D_ref","_prem_16"]

# labels=["PREM_3DC","S40RTS_3DC","S40RTS_1DC","GLAD_M25","PREM"]
# file_names=["_prem_3D","_3D","_1D","_glad","_synt"]

# labels=["S40RTS_3DC","S40RTS_3DC_QRFSI12","PREM_QL6","PREM"]
# file_names=["_3D","_S80_3D","_prem_16","_synt"]

labels=["PREM","S1500_3DC","S5000_3DC","S40RTS_3DC","S11000_3DC","S16000_3DC"]
file_names=["_synt","_3D_1","_3D_5","_3D","_3D_11","_3D_18"]

SMALL_SIZE = 10
MEDIUM_SIZE = 14
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labelsc
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)


handles=[]
for i,df in enumerate(file_names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker=markers[i], label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends")
plt.close()

phase_list=['S','SS','SSS','Sdiff']
ks        = [0,1,2,4]

fig,ax=plt.subplots(figsize=(7.08,4))
for j in range(len(phase_list)):
    ph=phase_list[j]
    k=ks[j]
    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all_real.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))]
    df_plot=df_plot[(df_plot["7"]>30)]
    #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(all_f):
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+df])<20)]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+1)+df])>0.7)]
        #df_plot=df_plot[df_plot[str(18+k*11+5)+df].notna()]
        #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+4)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+3)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+3)+df])<1.5)]
        #df_plot=df_plot[df_plot[str(18+k*11+5)+df].notna()]
        #df_plot=df_plot[(np.abs(np.log(df_plot[str(18+k*11+5)+df])<1.5))]
        df_plot=df_plot[df_plot[str(18+k*11+6)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+6)+df])<1.5)]
    df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    #df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]
    print(len(df_plot["1"]))
    # df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    # ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    # ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    # ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    # (df_plot["1"].str.split(".").str[0].isin(gsn))]
    # print(len(df_plot["1"]))



    if k==0:    
        df_plot = df_plot[((df_plot["7"]>30) & (df_plot["7"]<75)) | ((df_plot["7"]>85) & (df_plot["7"]<100))]
        ranges  = [30,60,100]
    elif k==1:
        df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
        ranges  = [50,80,110,140]
    elif k==2:
        df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<150))] #| ((df_plot["7"]>112) & (df_plot["7"]<150))]
        ranges  = [70,100,125,150]
    elif k==3:
        df_plot = df_plot[((df_plot["7"]>5) & (df_plot["7"]<25)) | ((df_plot["7"]>50) & (df_plot["7"]<65))]
        ranges=[5,65]
    elif k==4:
        df_plot = df_plot[(df_plot["7"]>100) & (df_plot["7"]<165)]
        ranges  = [100,125,150]
    
    misfits=[]
    i_s=[]
    for i,df in enumerate(file_names):
        sum_n=(df_plot[str(18+k*11+6)+df]**2).mean()
        mf   =0.5*sum_n
        misfits.append(mf)
        i_s.append(i)

    ax.plot(i_s[1:],misfits[1:],linestyle="--",marker=markers[k],color=colors[k],label=ph+" "+str(np.round(misfits[0],3)),lw=1)
ax.set_xticks(i_s[1:])
ax.set_xticklabels(labels[1:])
ax.set_xlabel("Global Models")
ax.set_ylabel("Average Misfit")
#ax.set_ylim(0.5,1.5)
#ax.set_title("Normalized with PREM Misfit="+str(misfits[0]))
ax.legend()
ax.grid()
plt.savefig(plot_dir+"/misfit_red_f_df.png")


fig,ax=plt.subplots(figsize=(7.08,4))
for j in range(len(phase_list)):
    ph=phase_list[j]
    k=ks[j]
    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all_real.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))]
    df_plot=df_plot[(df_plot["7"]>30)]
    #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(all_f):
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+df])<20)]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+1)+df])>0.7)]
        #df_plot=df_plot[df_plot[str(18+k*11+5)+df].notna()]
        #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+4)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+3)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+3)+df])<1.5)]
        # df_plot=df_plot[df_plot[str(18+k*11+5)+df].notna()]
        # df_plot=df_plot[(np.abs(np.log(df_plot[str(18+k*11+5)+df])<1.5))]
        df_plot=df_plot[df_plot[str(18+k*11+6)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+6)+df])<1.5)]
    df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    #df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]
    print(len(df_plot["1"]))
    # df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    # ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    # ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    # ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    # (df_plot["1"].str.split(".").str[0].isin(gsn))]
    # print(len(df_plot["1"]))



    if k==0:    
        df_plot = df_plot[((df_plot["7"]>30) & (df_plot["7"]<75)) | ((df_plot["7"]>85) & (df_plot["7"]<100))]
        ranges  = [30,60,100]
    elif k==1:
        df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
        ranges  = [50,80,110,140]
    elif k==2:
        df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<150))] #| ((df_plot["7"]>112) & (df_plot["7"]<150))]
        ranges  = [70,100,125,150]
    elif k==3:
        df_plot = df_plot[((df_plot["7"]>5) & (df_plot["7"]<25)) | ((df_plot["7"]>50) & (df_plot["7"]<65))]
        ranges=[5,65]
    elif k==4:
        df_plot = df_plot[(df_plot["7"]>100) & (df_plot["7"]<165)]
        ranges  = [100,125,150]
    
    misfits=[]
    i_s=[]
    for i,df in enumerate(file_names):
        sum_n=(df_plot[str(18+k*11+2)+df]**2).mean()
        mf   =0.5*sum_n
        misfits.append(mf)
        i_s.append(i)

    ax.plot(i_s[1:],misfits[1:],linestyle="--",marker=markers[k],color=colors[k],label=ph+" "+str(np.round(misfits[0],3)),lw=1)
ax.set_xticks(i_s[1:])
ax.set_xticklabels(labels[1:])
ax.set_xlabel("Global Models")
ax.set_ylabel("Average Misfit")
#ax.set_title("Normalized with PREM Misfit="+str(misfits[0]))

ax.legend()
ax.grid()
plt.savefig(plot_dir+"/misfit_cc_f_df.png")
