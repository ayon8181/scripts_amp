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
colors=['peru','black','forestgreen','magenta','mediumpurple','deepskyblue','darkorange','grey','crimson']     
colors=["red","black","royalblue","royalblue","darkmagenta","forestgreen","red"]
lss   =[1,2,1.5,1.25,1,1,1,1,1,1,1,1,1,1] 
style =['solid','solid',"solid","solid","solid","solid","solid","solid","solid"]
labels=["Real Data","PREMc","S8000","S11000","S5000","SPREMc","GLAD_M25","1D_REF","PREM_QL6.txt","PREM_QRFSI12","S40RTS_PREM_Crust"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None","None","None","None","None"] 
markers=["v",".","+","o","*","x","^"]
alp   =[1,1,1,1,1,1,1]
df=[]
#colors=["red","black"]
fill=["red",None]
plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/figures/attenuation/histograms_all/"



ev_list=[]
with open("./all.txt","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])






labels=["PREM_3DC","S40RTS_3DC","S40RTS_1DC","GLAD_M25","PREM"]
file_names=["_prem_3D","_3D","_1D","_glad","_synt"]

all_f=["_prem_3D","_3D","_1D","_glad","_synt"]

# labels=["S40RTS_3DC","S40RTS_3DC_QRFSI12","PREM","PREM_Q16","1D_REF"]
# file_names=["_3D","_S80_3D","_synt","_prem_16","_1D_ref"]

# # labels=["S40RTS_3DC","S40RTS_1DC","GLAD_M25","PREM_3DC","PREM"]
# # file_names=["_3D","_1D","_glad","_prem_3D","_synt"]
# all_f=["_3D","_S80_3D","_prem_16","_synt","_1D_ref"]
# # labels=["PREM_3DC","S40RTS_3DC","S40RTS_1DC","GLAD_M25","PREM"]
# # file_names=["_prem_3D","_3D","_1D","_glad","_synt"]

labels=["PREM_QL6","PREM","1D_REF"]
file_names=["_prem_16","_synt","_1D_ref"]

labels=["S40RTS_3DC","S40RTS_3DC_QRFSI12"]
file_names=["_3D","_S80_3D"]


SMALL_SIZE = 20
MEDIUM_SIZE = 24
BIGGER_SIZE = 50

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

phase_list=['S','SS','SSS','ScS','Sdiff']
def plot(k,ph):
    print("run")
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
    print(len(df_plot["1"]))



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
     

    for q in range(len(ranges)):
        if q==0:
           df_plots=df_plot
        else:
             df_plots=df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ]
        means=[]
        std=[]
        ## Plot Line Histograms
        plt.figure(1,figsize=(7.08,3.54))
        for i,df in enumerate(file_names):
            
            if k!=2:
               y, edges = np.histogram(np.log(df_plots[str(18+k*11+5)+df]), bins=15, range=(-1.45,1.45))
            else:
                y, edges = np.histogram(np.log(df_plots[str(18+k*11+5)+df]), bins=15, range=(-1.45,1.45))
            centers = 0.5 * (edges[1:] + edges[:-1])
            #plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)\
            plt.stairs(y/np.sum(y),edges,color=colors[i],fill=fill[i],label=labels[i])
            means.append(np.mean(np.log(df_plots[str(18+k*11+5)+df])))
            std.append(np.std(np.log(df_plots[str(18+k*11+5)+df])))
        plt.text(1.0,0.8*max(y/np.sum(y)),ph)
        if k!=2:
           plt.xlim(-1.45,1.45)
        else:
            plt.xlim(-1.45,1.45)
        plt.legend(fontsize=12)
        plt.xlabel("$\\chi_{amp_2}$")
        plt.ylabel("Normalized Counts")
        plt.tight_layout()
        plt.grid()
        plt.legend(fontsize=12)
        plt.savefig(plot_dir+"/line_amp_r_"+ph+"_"+str(q)+".png",dpi=600)
        plt.close()
    
        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" X2 mean "+"\n")
             for m,mea in enumerate(means):
                 txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" for range "+str(q)+"\n")
             txt.write(str(np.sum(y))+" "+"\n")

    
        means=[]
        std=[]
        # Plot Histograms for Cross Correlation Traveltime
        plt.figure(1,figsize=(7.08,3.54))
        for i,df in enumerate(file_names):
            if k!=2:
               y, edges = np.histogram(df_plots[str(18+k*11+2)+df], bins=20, range=(-20,20))
            else:
                y, edges = np.histogram(df_plots[str(18+k*11+2)+df], bins=20, range=(-20,20))
            centers = 0.5 * (edges[1:] + edges[:-1])
            #plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)
            plt.stairs(y/np.sum(y),edges,color=colors[i],fill=fill[i],label=labels[i])
            means.append(np.mean(df_plots[str(18+k*11+2)+df]))
            std.append(np.std(df_plots[str(18+k*11+2)+df]))
        plt.text(16,0.8*max(y/np.sum(y)),ph)
        plt.xlim(-20,20)
        #plt.legend(fontsize=8)(fontsize=12)
        plt.xlabel("Cross-correlation Travel Time Shift")
        plt.ylabel("Normalized Counts")
        plt.tight_layout()
        plt.grid()
        #plt.legend(fontsize=12)
        plt.savefig(plot_dir+"/line_cc_t_"+ph+str(q)+".png",dpi=600)
        plt.close()
        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" cc mean "+"\n")
             for m,mea in enumerate(means):
                 txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" for range "+str(q)+"\n")
             txt.write(str(np.sum(y))+" "+"\n")
    
        means=[]
        std=[]
        # Plot Histograms for Envelope Ratios
        plt.figure(1,figsize=(7.08,3.54))
        for i,df in enumerate(file_names):
            if k!=2:
               y, edges = np.histogram(df_plots[str(18+k*11+3)+df], bins=15, range=(-1.45,1.45))
            else:
                y, edges = np.histogram(df_plots[str(18+k*11+3)+df], bins=15, range=(-1.45,1.45))
            centers = 0.5 * (edges[1:] + edges[:-1])
            #plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)
            plt.stairs(y/np.sum(y),edges,color=colors[i],fill=fill[i],label=labels[i])
            means.append(np.mean(np.log(df_plots[str(18+k*11+3)+df])))
            std.append(np.std(np.log(df_plots[str(18+k*11+3)+df])))
        plt.text(1.0,0.8*max(y/np.sum(y)),ph)
        if k!=2:
           plt.xlim(-1.45,1.45)
        else:
            plt.xlim(-1.45,1.45)
        #plt.legend(fontsize=8)(fontsize=12)
        plt.xlabel("$\\chi_{env}$")
        plt.ylabel("Normalized Counts")
        plt.tight_layout()
        plt.grid()
        #plt.legend(fontsize=12)
        plt.savefig(plot_dir+"/line_envp_3d_1d_"+ph+str(q)+".png",dpi=600)
        plt.close()

        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" X_env mean "+"\n")
             for m,mea in enumerate(means):
                 txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" for range "+str(q)+"\n")
             txt.write(str(np.sum(y))+" "+"\n")
        means=[]
        std=[]
        plt.figure(1,figsize=(7.08,3.54))
        for i,df in enumerate(file_names):
            if k!=2:
               y, edges = np.histogram(df_plots[str(18+k*11+6)+df], bins=15, range=(-1.45,1.45))
            else:
                y, edges = np.histogram(df_plots[str(18+k*11+6)+df], bins=15, range=(-1.45,1.45))
            centers = 0.5 * (edges[1:] + edges[:-1])
            #plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)
            plt.stairs(y/np.sum(y),edges,color=colors[i],fill=fill[i],label=labels[i])
            means.append(np.mean(df_plots[str(18+k*11+6)+df]))
            std.append(np.std(df_plots[str(18+k*11+6)+df]))
        plt.text(1.0,0.8*max(y/np.sum(y)),ph)
        if k!=2:
           plt.xlim(-1.45,1.45)
        else:
            plt.xlim(-1.45,1.45)
        #plt.legend(fontsize=8)(fontsize=12)
        plt.xlabel("$\\chi_{amp}$")
        plt.ylabel("Normalized Counts")
        plt.tight_layout()
        plt.grid()
        #plt.legend(fontsize=10, position="top left")
        plt.savefig(plot_dir+"/line_amp_misfit_"+ph+str(q)+".png",dpi=600)
        plt.close()

        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" X_1 mean "+"\n")
             for m,mea in enumerate(means):
                 txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" for range "+str(q)+"\n")
             txt.write(str(np.sum(y))+" "+"\n")
    """
    ##Plot Scatters
    lab = ["$\\chi_{amp}$","$\\chi_{env}$","$\\chi_{amp_2}$",]
    nam = ["amp","env","bs_amp",]
    for b,g in enumerate(range(6,7)):
        
        fig,ax=plt.subplots(1,1,figsize=(7.08,3.54))
        
        ax.plot(np.linspace(-1.45,1.45,100), np.linspace(-1.45,1.5,100), linestyle='--', color='black')
        ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        if b==2:
           ax.scatter(df_plot[str(18+k*11+g)+"_prem_3D"],df_plot[str(18+k*11+g)+"_3D"],marker=".",alpha=0.6)
        else:
            ax.scatter(np.log(df_plot[str(18+k*11+g)+"_prem_3D"]),np.log(df_plot[str(18+k*11+g)+"_3D"]),marker=".",alpha=0.6)
        ax.set_xlabel(lab[b]+" for PREM_3DC")
        ax.set_ylabel(lab[b]+" for S40RTS_3DC")
        cc=df_plot[str(18+k*11+g)+"_prem_3D"].corr(df_plot[str(18+k*11+g)+"_3D"])
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        plt.tight_layout()
        plt.savefig(plot_dir+"/Scatter_"+nam[b]+"_"+ph+"_S8000.png",dpi=600)
        plt.close()
        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" Correlation between PREM_3DC and S40RTS_3DC is "+str(cc)+" for "+lab[b]+"\n")


        fig,ax=plt.subplots(1,1,figsize=(7.08,3.54))
        
        ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        if b==2:
           ax.scatter(df_plot[str(18+k*11+g)+"_prem_3D"],df_plot[str(18+k*11+g)+"_glad"],marker=".",alpha=0.6)
        else:
            ax.scatter(np.log(df_plot[str(18+k*11+g)+"_prem_3D"]),np.log(df_plot[str(18+k*11+g)+"_glad"]),marker=".",alpha=0.6)
        ax.set_xlabel(lab[b]+" for PREM_3DC")
        ax.set_ylabel(lab[b]+" for GLAD_M25")
        cc=df_plot[str(18+k*11+g)+"_prem_3D"].corr(df_plot[str(18+k*11+g)+"_glad"])
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        plt.tight_layout()
        plt.savefig(plot_dir+"/Scatter_"+nam[b]+"_"+ph+"_glad.png",dpi=600)
        plt.close()
        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" Correlation between PREM_3DC and GLAD_M25 is "+str(cc)+" for "+lab[b]+"\n")
    
        fig,ax=plt.subplots(1,1,figsize=(7.08,3.54))
        
        ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        if b==2:
           ax.scatter(df_plot[str(18+k*11+g)+"_3D"],df_plot[str(18+k*11+g)+"_glad"],marker=".",alpha=0.6)
        else:
            ax.scatter(np.log(df_plot[str(18+k*11+g)+"_3D"]),np.log(df_plot[str(18+k*11+g)+"_glad"]),marker=".",alpha=0.6)
        ax.set_xlabel(lab[b]+" for S40RTS_3DC")
        ax.set_ylabel(lab[b]+" for GLAD_M25")
        cc=df_plot[str(18+k*11+g)+"_3D"].corr(df_plot[str(18+k*11+g)+"_glad"])
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        plt.tight_layout()
        plt.savefig(plot_dir+"/Scatter_"+nam[b]+"_"+ph+"_S8000_glad.png",dpi=600)
        plt.close()
        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" Correlation between S40RTS_3DC and GLAD_M25 is "+str(cc)+" for "+lab[b]+"\n")

        fig,ax=plt.subplots(1,1,figsize=(7.08,3.54))
        
        ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        if b==2:
           ax.scatter(df_plot[str(18+k*11+g)+"_1D"],df_plot[str(18+k*11+g)+"_3D"],marker=".",alpha=0.6)
        else:
            ax.scatter(np.log(df_plot[str(18+k*11+g)+"_1D"]),np.log(df_plot[str(18+k*11+g)+"_3D"]),marker=".",alpha=0.6)
        ax.set_xlabel(lab[b]+" for S40RTS_1DC")
        ax.set_ylabel(lab[b]+" for S40RTS_3DC")
        cc=df_plot[str(18+k*11+g)+"_1D"].corr(df_plot[str(18+k*11+g)+"_3D"])
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        plt.tight_layout()
        plt.savefig(plot_dir+"/Scatter_"+nam[b]+"_"+ph+"_3D_1D.png",dpi=600)
        plt.close()
        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" Correlation between S40RTS_3DC and S40RTS_1DC is "+str(cc)+" for "+lab[b]+"\n")

        fig,ax=plt.subplots(1,1,figsize=(7.08,3.54))
        
        ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
        if b==2:
           ax.scatter(df_plot[str(18+k*11+g)+"_1D"],df_plot[str(18+k*11+g)+"_3D"],marker=".",alpha=0.6)
        else:
            ax.scatter(np.log(df_plot[str(18+k*11+g)+"_1D"]),np.log(df_plot[str(18+k*11+g)+"_3D"]),marker=".",alpha=0.6)
        ax.set_xlabel(lab[b]+" for S40RTS_1DC")
        ax.set_ylabel(lab[b]+" for PREM_3DC")
        cc=df_plot[str(18+k*11+g)+"_1D"].corr(df_plot[str(18+k*11+g)+"_prem_3D"])
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        plt.tight_layout()
        plt.savefig(plot_dir+"/Scatter_"+nam[b]+"_"+ph+"_prem_3D_1D.png",dpi=600)
        plt.close()
        with open(plot_dir+"/outputs.txt","a") as txt:
             txt.write(ph+" Correlation between PREM_3DC and S40RTS_1DC is "+str(cc)+" for "+lab[b]+"\n")
    """
    
for k,ph in enumerate(phase_list):
    plot(k,ph)
    

# pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
# results=pool.starmap(plot, [(k,ph)  for k,ph in enumerate(phase_list)])
# pool.close()
# pool.join()
