#!/usr/bin/env python
import numpy as np
import matplotlib.ticker as ticker
import csv
import pandas as pd
import os
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
import math
import multiprocessing as mp
import csv
from obspy.geodetics import locations2degrees
import scipy.stats as stats
from matplotlib.lines import Line2D

def calc_epdist(lon1, lat1, lon2, lat2):
    # Calculate the forward and back azimuths and distance
    dist = locations2degrees(lat1, lon1, lat2, lon2)
    return dist



def weight(event,df,t):
    
    df=df[df["0"]==event]
    #print(df)
    num_st = len(df)
    if num_st==0:
        return None
    weights=[]
    sta_list=df["1"].to_numpy()
    latitudes = df["6"].to_numpy()
    longitudes = df["5"].to_numpy()
    
   
    # Compute pairwise distances
    lat1, lat2 = np.meshgrid(latitudes, latitudes)
    lon1, lon2 = np.meshgrid(longitudes, longitudes)
    dist_matrix = calc_epdist(lon1, lat1, lon2, lat2)  # Ensure calc_epdist supports arrays
    
    # Apply the exponential decay function element-wise
    A = np.exp(-(dist_matrix / t) ** 2)
        #print(A)
        
        
                
    weights.append(A.sum(axis=1))

    # print(weights)
    # print(sta_list)
    return [event,sta_list,weights,t]

us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"] 
file_names      = ["_real","_prem_3D","_3D","_3D_2","_3D_15","_1D","_glad","_1D_ref","_prem_q16","_3D_atten","_1D"]
colors=["dodgerblue","red","black", 'blue','mediumpurple','deepskyblue','darkorange','grey','crimson']   
lss   =[2,1.5,1.75,1.25,0.5,1,1,1,1,1,1,1,1,1] 
style =['solid','solid',"solid","solid","solid","solid","solid","solid","solid"]
labels=["Real Data","PREMc","S8000","S11000","S5000","SPREMc","GLAD_M25","1D_REF","PREM_QL6.txt","PREM_QRFSI12","S40RTS_PREM_Crust"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None","None","None","None","None"] 
markers=["v",".","+","o","*","x","^"]
alp   =[1,1,1,1,2,1,1]
df=[]

plot_dir="/home/ayon/WORK/grl/outputs2/new_figures/relative/"



ev_list=[]
with open("event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])









#all_f=["_glad","_3D","_S80_3D","_synt","_1D_ref","_prem_16"]

# labels=["PREM_3DC","S40RTS_3DC","S40RTS_1DC","GLAD_M25","PREM"]
# file_names=["_prem_3D","_3D","_1D","_glad","_synt"]

# labels=["S40RTS_3DC","S40RTS_3DC_QRFSI12","PREM_QL6","PREM"]
# file_names=["_3D","_S80_3D","_prem_16","_synt"]
file_names=[]
labels=[]

n_list=[872,1482,2405,3326,3941,4534,4922,5405,6026,6868,8112,9010,11289,12181,13332,15821,15988]
for n in n_list:
    labels.append(str(n))
    file_names.append("_"+str(n))


file_names.append("_mod_3")
file_names.append("_glad")
file_names.append("_real")
#file_names2=file_names[0:-2]


SMALL_SIZE = 28
MEDIUM_SIZE = 28
BIGGER_SIZE = 50

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labelsc
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

font = {'family' : 'serif',
         'serif':  'cmr10'
         }

plt.rc('font', **font,size=SMALL_SIZE)  
plt.rc('axes.formatter', use_mathtext=True)
plt.rc("mathtext", fontset="cm")


handles=[]
phase_list=['55-110','110-140']
for i,df in enumerate(phase_list):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker="o", label=df))
                                    
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends2")
plt.close()

phase_list=['SS_S','SS_Sdiff']
ks        = [0,1,2,4]
total_mf_amp=np.zeros(len(n_list))
total_mf_cc=np.zeros(len(n_list))
total_mf_wf=np.zeros(len(n_list))


fig,(ax,ax2)=plt.subplots(1,2,figsize=(16,7.5))
fig.text(0.029, 0.9, 'a)', fontsize=28, fontweight='bold')  # Position of "a)"
fig.text(0.5, 0.9, 'b)', fontsize=28, fontweight='bold') 

fig_hist,ax_hist=plt.subplots(2,2,figsize=(14,14))
fig_hist.text(0.029, 0.9, 'a)', fontsize=28, fontweight='bold')  # Position of "a)"
fig_hist.text(0.5, 0.9, 'b)', fontsize=28, fontweight='bold')
fig_hist.text(0.029, 0.5, 'c)', fontsize=28, fontweight='bold')  # Position of "a)"
fig_hist.text(0.5, 0.5, 'd)', fontsize=28, fontweight='bold')

ax_hist[0,0].set_ylabel("Counts")
ax_hist[1,0].set_ylabel("Counts")
ax_hist[0,0].set_xlabel(r"$\delta \rm lnA$")
ax_hist[0,1].set_xlabel(r"$\delta \rm lnA$")
ax_hist[1,0].set_xlabel(r"$\delta \rm lnA$ residual")
ax_hist[1,1].set_xlabel(r"$\delta \rm lnA$ residual")
    


df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
#df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]

df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<110)]
for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]


misfits=[]
means=[]
stds=[]
i_s=[]
counts=[]
len_n=0
total_glad=0
df_plot["dlnA"+"_real"]=df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"]
df_plot["dlnA"+"_glad"]=df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"]
#nrm=((df_plot["dlnA"+"_real"]).std())**2*len(df_plot)
nrm=np.sum(df_plot["dlnA"+"_real"]**2)
for i,df in enumerate(file_names[:-3]):
    df_plot["dlnA"+df]=df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df]
    #sum_n=((1.0/(df_plot['weights']*df_plot['fow']))*(df_plot[str(8)+df])**2).sum()
    #res_v=((df_plot["dlnA"+df]-df_plot["dlnA"+"_real"]).std())**2
    res_v=np.sum((df_plot["dlnA"+df]-df_plot["dlnA"+"_real"])**2)
    len_n=len(df_plot)
    misfits.append(res_v)
    i_s.append(n_list[i])

y,edges=np.histogram(df_plot["dlnA"+"_8112"],bins=15,range=(-2,2))
ax_hist[0,0].stairs(y/np.sum(y), edges, color="red",label="S40RTS",fill=False,lw=2)
mean,std = np.mean(df_plot["dlnA"+"_8112"]),np.std(df_plot["dlnA"+"_8112"])
# ax_hist[0,0].axvline(mean, color='red', linestyle='dashed', linewidth=1)
# ax_hist[0,0].axvline(mean+std, color='red', linestyle='dotted', linewidth=1)
# ax_hist[0,0].axvline(mean-std, color='red', linestyle='dotted', linewidth=1)
ax_hist[0,0].text(0.05,0.9,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[0,0].transAxes,color="red",fontsize=16)
y,edges=np.histogram(df_plot["dlnA"+"_glad"],bins=15,range=(-2,2))
ax_hist[0,0].stairs(y/np.sum(y), edges, color="dodgerblue",label="GLAD-M25",fill=False,lw=2)
mean,std = np.mean(df_plot["dlnA"+"_glad"]),np.std(df_plot["dlnA"+"_glad"])
# ax_hist[0,0].axvline(mean, color='dodgerblue', linestyle='dashed', linewidth=1)
# ax_hist[0,0].axvline(mean+std, color='dodgerblue', linestyle='dotted', linewidth=1)
# ax_hist[0,0].axvline(mean-std, color='dodgerblue', linestyle='dotted', linewidth=1)
ax_hist[0,0].text(0.05,0.8,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[0,0].transAxes,color="dodgerblue",fontsize=16)
y,edges=np.histogram(df_plot["dlnA"+"_real"],bins=15,range=(-2,2))
ax_hist[0,0].stairs(y/np.sum(y), edges, color="black",label="Recorded",fill=False,lw=2)
mean,std = np.mean(df_plot["dlnA"+"_real"]),np.std(df_plot["dlnA"+"_real"])
# ax_hist[0,0].axvline(mean, color='black', linestyle='dashed', linewidth=1)    
# ax_hist[0,0].axvline(mean+std, color='black', linestyle='dotted', linewidth=1)
# ax_hist[0,0].axvline(mean-std, color='black', linestyle='dotted', linewidth=1)
ax_hist[0,0].text(0.05,0.7,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[0,0].transAxes,color="black",fontsize=16)
ax_hist[0,0].grid()
y,edges=np.histogram(-(df_plot["dlnA"+"_8112"]-df_plot["dlnA"+"_real"]),bins=15,range=(-2,2))
ax_hist[1,0].stairs(y/np.sum(y), edges, color="red",fill=False,lw=2)
mean,std = np.mean(-(df_plot["dlnA"+"_8112"]-df_plot["dlnA"+"_real"])),np.std(-(df_plot["dlnA"+"_8112"]-df_plot["dlnA"+"_real"]))
# ax_hist[1,0].axvline(mean, color='red', linestyle='dashed', linewidth=1)
# ax_hist[1,0].axvline(mean+std, color='red', linestyle='dotted', linewidth=1)
# ax_hist[1,0].axvline(mean-std, color='red', linestyle='dotted', linewidth=1)
ax_hist[1,0].text(0.05,0.9,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[1,0].transAxes,color="red",fontsize=16)
y,edges=np.histogram(-(df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]),bins=15,range=(-2,2))
ax_hist[1,0].stairs(y/np.sum(y), edges, color="dodgerblue",fill=False,lw=2)
mean,std = np.mean(-(df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"])),np.std(-(df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]))
# ax_hist[1,0].axvline(mean, color='dodgerblue', linestyle='dashed', linewidth=1)
# ax_hist[1,0].axvline(mean+std, color='dodgerblue', linestyle='dotted', linewidth=1)
# ax_hist[1,0].axvline(mean-std, color='dodgerblue', linestyle='dotted', linewidth=1)
ax_hist[1,0].text(0.05,0.8,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[1,0].transAxes,color="dodgerblue",fontsize=16)
ax_hist[1,0].grid()



#res_glad=((df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]).std())**2 
res_glad=np.sum((df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"])**2)



ax2.plot(i_s,misfits/nrm,linestyle="--",marker='o',color=colors[0],lw=2,label=r"$55-110\degree$")
ax2.scatter(8112,res_glad/nrm,marker="*",color=colors[0],s=200)
#ax2.plot(0,sum_real,linestyle="--",marker="^",color=colors[0])
print(len(df_plot))
    #ax2.plot(i_s[1:],np.array(counts)/len(df_plot),linestyle="--",marker=markers[j],color=colors[j],label=ph,lw=1)

df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]
#df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]

df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]
for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]


misfits=[]
means=[]
stds=[]
i_s=[]
counts=[]
print(len(df_plot))
df_plot["dlnA"+"_real"]=df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"]
df_plot["dlnA"+"_glad"]=df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"]
#nrm=((df_plot["dlnA"+"_real"]).std())**2*len(df_plot)
nrm=np.sum(df_plot["dlnA"+"_real"]**2)
for i,df in enumerate(file_names[:-3]):
    df_plot["dlnA"+df]=df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df]
    #sum_n=((1.0/(df_plot['weights']*df_plot['fow']))*(df_plot[str(8)+df])**2).sum()
    res_v=np.sum((df_plot["dlnA"+df]-df_plot["dlnA"+"_real"])**2)
    #res_v=((df_plot["dlnA"+df]-df_plot["dlnA"+"_real"]).std())**2
    len_n=len(df_plot)
    misfits.append(res_v)
    i_s.append(n_list[i])


y,edges=np.histogram(df_plot["dlnA"+"_8112"],bins=15,range=(-2,2))
ax_hist[0,1].stairs(y/np.sum(y), edges, color="red",label="S40RTS",fill=False,lw=2)
mean,std=np.mean(df_plot["dlnA"+"_8112"]),np.std(df_plot["dlnA"+"_8112"])
# ax_hist[0,1].axvline(mean, color='red', linestyle='dashed', linewidth=1)
# ax_hist[0,1].axvline(mean+std, color='red', linestyle='dotted', linewidth=1)
# ax_hist[0,1].axvline(mean-std, color='red', linestyle='dotted', linewidth=1)
ax_hist[0,1].text(0.05,0.9,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[0,1].transAxes,color="red",fontsize=16)
y,edges=np.histogram(df_plot["dlnA"+"_glad"],bins=15,range=(-2,2))
ax_hist[0,1].stairs(y/np.sum(y), edges, color="dodgerblue",label="GLAD-M25",fill=False,lw=2)
mean,std=np.mean(df_plot["dlnA"+"_glad"]),np.std(df_plot["dlnA"+"_glad"])
# ax_hist[0,1].axvline(mean, color='dodgerblue', linestyle='dashed', linewidth=1)
# ax_hist[0,1].axvline(mean+std, color='dodgerblue', linestyle='dotted', linewidth=1)
# ax_hist[0,1].axvlie(mean-std, color='dodgerblue', linestyle='dotted', linewidth=1)
ax_hist[0,1].text(0.05,0.8,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[0,1].transAxes,color="dodgerblue",fontsize=16)
y,edges=np.histogram(df_plot["dlnA"+"_real"],bins=15,range=(-2,2))
ax_hist[0,1].stairs(y/np.sum(y), edges, color="black",label="Recorded",fill=False,lw=2)
mean,std=np.mean(df_plot["dlnA"+"_real"]),np.std(df_plot["dlnA"+"_real"])
ax_hist[0,1].text(0.05,0.7,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[0,1].transAxes,color="black",fontsize=16)
ax_hist[0,1].grid()
y,edges=np.histogram(-(df_plot["dlnA"+"_8112"]-df_plot["dlnA"+"_real"]),bins=15,range=(-2,2))
ax_hist[1,1].stairs(y/np.sum(y), edges, color="red",fill=False,lw=2)
mean,std=np.mean(-(df_plot["dlnA"+"_8112"]-df_plot["dlnA"+"_real"])),np.std(-(df_plot["dlnA"+"_8112"]-df_plot["dlnA"+"_real"]))
# ax_hist[1,1].axvline(mean, color='red', linestyle='dashed', linewidth=1)
# ax_hist[1,1].axvline(mean+std, color='red', linestyle='dotted', linewidth=1)
# ax_hist[1,1].axvline(mean-std, color='red', linestyle='dotted', linewidth=1)
ax_hist[1,1].text(0.05,0.9,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[1,1].transAxes,color="red",fontsize=16)
y,edges=np.histogram(-(df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]),bins=15,range=(-2,2))
ax_hist[1,1].stairs(y/np.sum(y), edges, color="dodgerblue",fill=False,lw=2)
mean,std=np.mean(-(df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"])),np.std(-(df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]))
# ax_hist[1,1].axvline(mean, color='dodgerblue', linestyle='dashed', linewidth=1)
# ax_hist[1,1].axvline(mean+std, color='dodgerblue', linestyle='dotted', linewidth=1)
# ax_hist[1,1].axvline(mean-std, color='dodgerblue', linestyle='dotted', linewidth=1)
ax_hist[1,1].text(0.05,0.8,"Mean: "+str(round(mean,2))+" Std: "+str(round(std,2)),transform=ax_hist[1,1].transAxes,color="dodgerblue",fontsize=16)
ax_hist[1,1].grid()
ax_hist[0,0].set_xlim(-2,2)
ax_hist[0,1].set_xlim(-2,2)
ax_hist[1,0].set_xlim(-2,2)
ax_hist[1,1].set_xlim(-2,2)
fig_hist.tight_layout()
ax_hist[0,1].legend(fontsize=16)
fig_hist.savefig(plot_dir+"/"+"hist.pdf",dpi=600)


#res_glad=(((df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]).std())**2)
res_glad=np.sum((df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"])**2)
ax2.plot(i_s,misfits/nrm,linestyle="--",marker='o',color=colors[1],lw=2,label=r"$110-140\degree$")
ax2.scatter(8112,res_glad/nrm,marker="*",color=colors[1],s=200)







df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]
#df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]

df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<140)]
for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]

misfits=[]
means=[]
stds=[]
i_s=[]
counts=[]
print(len(df_plot))
df_plot["dlnA"+"_real"]=df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"]
df_plot["dlnA"+"_glad"]=df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"]
# nrm=((df_plot["dlnA"+"_real"]).std())**2*len(df_plot)
nrm=np.sum(df_plot["dlnA"+"_real"]**2)
for i,df in enumerate(file_names[:-3]):
    df_plot["dlnA"+df]=df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df]
    #sum_n=((1.0/(df_plot['weights']*df_plot['fow']))*(df_plot[str(8)+df])**2).sum()
    #res_v=(((df_plot["dlnA"+df]-df_plot["dlnA"+"_real"]).std())**2)
    res_v=np.sum((df_plot["dlnA"+df]-df_plot["dlnA"+"_real"])**2)
    len_n=len(df_plot)
    misfits.append(res_v)
    i_s.append(n_list[i])

#res_glad=(((df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]).std())**2)
res_glad=np.sum((df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"])**2) 
ax2.plot(i_s,misfits/nrm,marker='o',color="black",lw=2,label=r"All")
ax2.scatter(8112,res_glad/nrm,marker="*",color="black",s=200)

print(i_s,misfits/nrm)
#ax.set_xticklabels(labels)
ax2.set_xlabel("N", fontsize=32)
ax2.set_ylabel(r"$\rm K^2_{ss/s}$",fontsize=32)
#ax.set_ylim(0.8,1.1)
#ax.set_ylim(0.8,1.2)
ax2.set_title("Amplitudes")
#ax.legend()
ax2.grid()
#plt.tight_layou
# t()
#fig.savefig(plot_dir+"/"+"20_40_misfit_relative_body.pdf",dpi=600)
ax2.set_xticks(np.arange(0, 16001, 4000))
ax2.set_xlim(0,16000)
#ax.set_yticks(np.arange(0.875, 1.076, 0.1))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(2000))


bin_edges=np.arange(55,146,10)
bin_centers = bin_edges[:-1] + 5
            
df_plot["bin"] = np.digitize(df_plot["7"], bin_edges)
df_plot["diff1"] = df_plot["dlnA"+"_glad"]-df_plot["dlnA"+"_real"]
df_plot["diff2"] = df_plot["dlnA"+"_8112"]-df_plot["dlnA"+"_real"]


# Compute mean and std deviation per bin
mean_diff1 = df_plot.groupby("bin")["diff1"].mean()
std_diff1 = df_plot.groupby("bin")["diff1"].std()

mean_diff2 = df_plot.groupby("bin")["diff2"].mean()
std_diff2 = df_plot.groupby("bin")["diff2"].std()

print(len(mean_diff1), bin_centers)
fig_dist,ax_dist=plt.subplots(1,1,figsize=(14,8))
bin_centers = bin_edges[:-1] + 5

ax_dist.errorbar(bin_centers, mean_diff1, yerr=std_diff1, fmt="o", color="dodgerblue", label="GLAD-M25", capsize=3,capthick=1,elinewidth=1,ms=16)
ax_dist.errorbar(bin_centers-1, mean_diff2, yerr=std_diff2, fmt="o", color="red", label="S40RTS", capsize=3,capthick=1,elinewidth=1,ms=16)
ax_dist.set_xlabel("Distance (degrees)", fontsize=32)
ax_dist.set_ylabel(r"$\delta \rm lnA$ residual", fontsize=32)
ax_dist.legend()
ax_dist.grid()
fig_dist.tight_layout()
fig_dist.savefig(plot_dir+"/"+"dist.pdf",dpi=600)
plt.close(fig_dist)

#fig.savefig("test.png",dpi=600)

#fig,ax=plt.subplots(figsize=(7.08,4))

df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]
#df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]

df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<110)]
for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]

misfits=[]
means=[]
stds=[]
i_s=[]
counts=[]
total_glad=0
len_n=0
nrm_total=0
df_plot["dT"+"_real"]=df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"]
df_plot["dT"+"_glad"]=df_plot[str(14+1*5+2)+"_glad"]-df_plot[str(14+0*5+2)+"_glad"]
nrm=np.mean(df_plot["dT"+"_real"]**2)
for i,df in enumerate(file_names[:-3]):
    df_plot["dT"+df]=df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df]
    #sum_n=((1.0/(df_plot['weights']*df_plot['fow']))*(df_plot[str(8)+df])**2).sum()
    #res_v=(((df_plot["dT"+df]-df_plot["dT"+"_real"]).std())**2)
    res_v=np.mean((df_plot["dT"+df]-df_plot["dT"+"_real"])**2)
    #print(misfits)
    misfits.append(res_v)
    i_s.append(n_list[i])
#sum_real=((df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"])**2).mean()
#res_glad=(((df_plot["dT"+"_glad"]-df_plot["dT"+"_real"]).std())**2) 
res_glad=np.mean((df_plot["dT"+"_glad"]-df_plot["dT"+"_real"])**2)  
print(len(df_plot))
ax.plot(i_s,misfits/nrm,linestyle="--",marker='o',color=colors[0],lw=2,label=r"$55-110\degree$")
ax.scatter(8112,res_glad/nrm,marker="*",color=colors[0],s=200)
    #ax2.plot(i_s[1:],np.array(counts)/len(df_plot),linestyle="--",marker=markers[j],color=colors[j],label=ph,lw=1)


df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]
#df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]

df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]

df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]
for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]

misfits=[]
means=[]
stds=[]
i_s=[]
counts=[]
df_plot["dT"+"_real"]=df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"]
df_plot["dT"+"_glad"]=df_plot[str(14+1*5+2)+"_glad"]-df_plot[str(14+0*5+2)+"_glad"]
#nrm=((df_plot["dT"+"_real"]).std())**2
nrm=np.mean(df_plot["dT"+"_real"]**2)
for i,df in enumerate(file_names[:-3]):
    df_plot["dT"+df]=df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df]
    #sum_n=((1.0/(df_plot['weights']*df_plot['fow']))*(df_plot[str(8)+df])**2).sum()
    #res_v=(((df_plot["dT"+df]-df_plot["dT"+"_real"]).std())**2)
    res_v=np.mean((df_plot["dT"+df]-df_plot["dT"+"_real"])**2)
    #print(misfits)
    misfits.append(res_v)
    i_s.append(n_list[i])
#sum_real=((df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"])**2).mean()
#res_glad=(((df_plot["dT"+"_glad"]-df_plot["dT"+"_real"]).std())**2) 
res_glad=np.mean((df_plot["dT"+"_glad"]-df_plot["dT"+"_real"])**2)
#print(misfits/nrm,i_s)
ax.plot(i_s,misfits/nrm,linestyle="--",marker='o',color=colors[1],lw=2,label=r"$110-140\degree$")
ax.scatter(8112,res_glad/nrm,marker="*",color=colors[1],s=200)


df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
df_plot=df_plot[(df_plot["0"].isin(ev_list))]

df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<140)]
for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]

misfits=[]
means=[]
stds=[]
i_s=[]
counts=[]
df_plot["dT"+"_real"]=df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"]
df_plot["dT"+"_glad"]=df_plot[str(14+1*5+2)+"_glad"]-df_plot[str(14+0*5+2)+"_glad"]
nrm=np.mean(df_plot["dT"+"_real"]**2)
for i,df in enumerate(file_names[:-3]):
    df_plot["dT"+df]=df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df]
    #sum_n=((1.0/(df_plot['weights']*df_plot['fow']))*(df_plot[str(8)+df])**2).sum()
    #res_v=(((df_plot["dT"+df]-df_plot["dT"+"_real"]).std())**2)
    res_v=np.mean((df_plot["dT"+df]-df_plot["dT"+"_real"])**2)
    #print(misfits)
    misfits.append(res_v)
    i_s.append(n_list[i])
#sum_real=((df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"])**2).mean()
#res_glad=(((df_plot["dT"+"_glad"]-df_plot["dT"+"_real"]).std())**2) 
res_glad=np.mean((df_plot["dT"+"_glad"]-df_plot["dT"+"_real"])**2)

ax.plot(i_s,misfits/nrm,marker='o',color="black",lw=2,label=r"All")
ax.scatter(8112,res_glad/nrm,marker="*",color="black",s=200)

print(i_s,misfits/nrm)
# df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
# #df_plot=df[df[str(14+k*5+5)+"_real"].notna()]
# df_plot=df[(df["0"].isin(ev_list))]
# df_plot=df_plot[(df_plot["7"]>=100) & (df_plot["7"]<140)]
# #df_plot=df_plot[(np.abs(df_plot[str(14+k*5+2)+"_real"])<15) & (np.abs(df_plot[str(14+k*5+2)+"_1D"])<15) &(np.abs(df_plot[str(14+k*5+2)+"_3D"])<15) & (np.abs(df_plot[str(14+k*5+2)+"_prem_3D"])<15)]
# for i,df in enumerate(file_names):
#     for k in [0,1]:
#         df_plot=df_plot[(np.abs(df_plot[str(14+k*5+2)+df])<12)]
#         df_plot=df_plot[(np.abs(df_plot[str(14+k*5+1)+df])>0.8)]
#         df_plot=df_plot[((np.log(np.abs(df_plot[str(14+k*5+4)+df])))<2)]

# #df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
# df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
#         ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
#         ~(df_plot["1"].str.split(".").str[0].isin(us)) &
#         ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
#         (df_plot["1"].str.split(".").str[0].isin(gsn))]

# misfits=[]
# means=[]
# stds=[]
# i_s=[]
# counts=[]
# for i,df in enumerate(file_names2[1:]):
#     #sum_n=((1.0/(df_plot['weights']*df_plot['fow']))*(df_plot[str(10)+df])**2).sum()
#     sum_n=((df_plot[str(8)+df])**2).sum()
#     sum_n=((df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df])**2).sum()
#     sum_prem=0.5*((df_plot[str(14+1*5+2)+"_prem"]-df_plot[str(14+0*5+2)+"_prem"])**2).sum()
#     sum_glad=0.5*((df_plot[str(14+1*5+2)+"_glad"]-df_plot[str(14+0*5+2)+"_glad"])**2).sum()
#     misfits.append(0.5*sum_n)
#     total_mf_amp[i]+=0.5*sum_n
#     i_s.append(n_list[i])

# ax.plot(i_s,misfits/misfits[0]linestyle="--",marker='o',color=colors[2],lw=1)
# ax.plot(8112,sum_glad,linestyle="--",marker="*",color=colors[2])
# ax.plot(8112,sum_real,linestyle="--",marker="*",color=colors[2])

#ax.set_xticklabels(labels)
ax.set_xlabel("N",fontsize=32)
ax.set_ylabel(r"$\rm K^2_{ss-s}$", fontsize=32)
#ax.set_ylim(0.8,1.1)
#ax.set_ylim(0.8,1.2)
ax.set_title("Traveltimes")
ax.set_xticks(np.arange(0, 16001, 4000))
ax.set_xlim(0,16000)
ax.set_ylim(0.2,0.8)
ax.set_yticks(np.arange(0.2, 0.8, 0.2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(2000))
#ax.legend()
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax2.set_yticks(np.arange(0.6, 1.1, 0.1))
ax2.set_ylim(0.6,1.1)

ax.grid(True, which='both', linestyle='dotted')
ax2.grid(True, which='both', linestyle='dotted')
ax2.legend()
plt.tight_layout()
fig.savefig(plot_dir+"/"+"residual_misfit.pdf",dpi=600)

