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
from importlib import reload


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
with open("./all.txt","r") as txt:
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


labels=["S40RTS_3DC","GLAD_M25","PREM","S40RTS_1DC"]#,"PREM_QL6","PREM_QRFSI12","S40RTS_3DC_QRFSI12"]
file_names=["_3D","_glad","_synt","_1D"]#,"_prem_16","_3D_atten","_S80_3D"]

phase_list=['S','SS','SSS','ScS','Sdiff']
def plot(k,ph):
    
    g=Geodesic.WGS84

    
    plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/figures/events/all"
    if os.path.exists(plot_dir)==False:
        os.mkdir(plot_dir)

    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all_real.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    #df_plot=df[(df["0"].isin(event))]
    df_plot=df[(df["7"]>30)]
    #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(file_names):
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+df])<18)]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+1)+df])<0.85)]
        #df_plot=df_plot[df_plot[str(18+k*11+3)+df].notna()]
        #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+3)+df])<1.5)]
        
        
        df_plot=df_plot[df_plot[str(18+k*11+6)+df].notna()]
        #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+6)+df])<1.5)]
    
    for i,df in enumerate(file_names):
        if k==0:    
            df_plots = df_plot[((df_plot["7"]>30) & (df_plot["7"]<75)) | ((df_plot["7"]>85) & (df_plot["7"]<100))]
        elif k==1:
            df_plots = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
        elif k==2:
            df_plots = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>115) & (df_plot["7"]<150))]
        elif k==3:
            df_plots = df_plot[((df_plot["7"]>5) & (df_plot["7"]<25)) | ((df_plot["7"]>50) & (df_plot["7"]<65))]
        elif k==4:
            df_plots = df_plot[(df_plot["7"]>100) & (df_plot["7"]<150)]
        df_plots=df_plots[df_plots[str(18+k*11+6)+df]<-0.8]
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='no_green',reverse=True,series=[-1.4,1.4])
        fig.basemap(region="d",projection="N-60/12c",frame=True)
        for v, row in df_plots.iterrows():
            #az=g.Inverse(row["3"],row["2"],row["6"],row["7"])['azi1']
            #print(az)
            fig.plot(x=[row["2"],row["5"]],y=[row["3"],row["6"]],zvalue=row[str(18+k*11+6)+df],pen="0.5p,+z",cmap=True)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar(frame=['a'+str(0.25)+'f0.25',f"x+l\"ln(A_3D/A_1D)\""])
        fig.savefig(plot_dir+"/amp_"+file_names[i]+"_"+ph+".png",dpi=600)
        plt.close()
    df_sorted=df_plots.sort_values(by=str(18+k*11+6)+df)[["0","1","2","3","4","5","6","7",str(18+k*11+6)+file_names[0],str(18+k*11+6)+file_names[1],str(18+k*11+6)+file_names[2],str(18+k*11+1)+file_names[0],str(18+k*11+1)+file_names[1],str(18+k*11+1)+file_names[2]]]
    df_sorted.to_csv(plot_dir+"/amp_file_"+"_"+ph+".csv")
    """
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='no_green',reverse=True,series=[-1.4,1.4])
        fig.basemap(region="d",projection="N-60/12c",frame=True)
        for v, row in df_plots.iterrows():
            #az=g.Inverse(row["3"],row["2"],row["6"],row["7"])['azi1']
            fig.plot(x=[row["2"],row["5"]],y=[row["3"],row["6"]],zvalue=row[str(18+k*11+3)+df],pen="0.5p,+z",cmap=True)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar(frame=['a'+str(0.25)+'f0.25',f"x+l\"ln(A_3D/A_1D)\""])
        fig.savefig(plot_dir+"/env_"+file_names[i]+"_"+ph+".png",dpi=600)
        plt.close()
        df_sorted=df_plots.sort_values(by=str(18+k*11+3)+df)[["0","1","2","3","4","5","6","7",str(18+k*11+3)+df]]
        df_sorted.to_csv(plot_dir+"/env_file_"+file_names[i]+"_"+ph+".csv")
     """    
#pygmt.exit()
"""
def plot_relative(event,k,ph):
    import pygmt
    reload(pygmt)
    g=Geodesic.WGS84

    
    plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/figures/events/all"
    if os.path.exists(plot_dir)==False:
        os.mkdir(plot_dir)

    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(event))]
    df_plot=df_plot[(df_plot["7"]>30)]
    if k==0:
        l=1
        m=0
        df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
        p=8
    elif k==1:
         l=1
         m=4
         df_plot = df_plot[(df_plot["7"]>100) & (df_plot["7"]<150)]
         p=10
    elif k==2:
         l=2
         m=1
         df_plot = df_plot[(((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<140)))]
         p=12
    elif k==3:
         l=2
         m=4
         df_plot = df_plot[(((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<150)))]
         p=16
    df_plots=df_plot
    for i,df in enumerate(file_names):
        df_plots=df_plots[(np.abs(df_plots[str(18+l*11+2)+df])<20)]
        df_plots=df_plots[(np.abs(df_plots[str(18+l*11+1)+df])>0.7)]
        df_plots=df_plots[(np.abs(df_plots[str(18+m*11+2)+df])<20)]
        df_plots=df_plots[(np.abs(df_plots[str(18+m*11+1)+df])>0.7)]
        print(len(df_plots))
        df_plots=df_plots[df_plots[str(p)+df].notna()]
        df_plots=df_plots[(np.abs(np.log((df_plots[str(p)+df])))<1.5)]
        print(len(df_plots))
    
    for i,df in enumerate(file_names):
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='no_green',reverse=True,series=[-1.4,1.4])
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        for v, row in df_plots.iterrows():
            #az=g.Inverse(row["3"],row["2"],row["6"],row["7"])['azi1']
            #print(az)
            fig.plot(x=[row["2"],row["5"]],y=[row["3"],row["6"]],zvalue=np.log(row[str(p)+df]),pen="0.5p,+z",cmap=True)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar(frame=['a'+str(0.25)+'f0.25',f"x+l\"ln("+ph+")\""])
        fig.savefig(plot_dir+"/amp_"+file_names[i]+"_"+ph+".png",dpi=600)
        plt.close()
        df_sorted=df_plots.sort_values(by=str(p)+df)[["0","1","2","3","4","5","6","7",str(p)+df]]
        df_sorted.to_csv(plot_dir+"/amp_file_"+file_names[i]+"_"+ph+".csv")

for k,ph in enumerate(phase_list):
    pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
    results=pool.starmap(plot, [(ev,k,ph)  for ev in ev_list])
    pool.close()
    pool.join()
ratio_list=["SS_S","SS_Sdiff","SSS_SS","SSS_Sdiff"]
#for k,ph in enumerate(ratio_list):

pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot, [(ev_list,k,ph)  for k,ph in enumerate(phase_list)])
pool.close()
pool.join()

pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot_relative, [(ev_list,k,ph)  for k,ph in enumerate(ratio_list)])
pool.close()
pool.join()
"""

plot(4,'Sdiff')
plot(0,'S')

    
    
        


