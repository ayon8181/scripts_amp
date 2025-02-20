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
import matplotlib.gridspec as gridspec
#import seaborn as sns
import scipy.stats as stats
import geopy.distance
from geographiclib.geodesic import Geodesic
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
#from mpi4py import MPI

font = { "family":"serif",
          "color": "darkred",
          "weight":"normal",
          }

SMALL_SIZE = 28
MEDIUM_SIZE = 32
BIGGER_SIZE = 50

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

font = {'family' : 'serif',
         'serif':  'cmr10'
         }

plt.rc('font', **font,size=SMALL_SIZE)  
plt.rc('axes.formatter', use_mathtext=True)


plt.rcParams['mathtext.default'] = 'rm'
plt.rc("mathtext", fontset="cm")
ev_list=[]
with open("event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])




us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"]         

file_names      = ["_real","_1D","_3D","_3D_2","_3D_15","_prem_3D","_glad","_synt",]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['olive','crimson','dodgerblue','darkgreen','black','peru','indianred','grey','black']    
lss   =[1,1,1,1,1,1,1,1]
style =['None','solid',"solid","solid","solid","solid",'solid']
labels=["real_data",'SPREMc','S8000',"S11000","S5000","PREMc","GLAD_M25"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None"] 
markers=[".","+","o","*","x","^","v","8",">","<","6"]
alp   =[0.5,1,1,1,1,1,1]
df=[]
plot_dir="new_figures_test"
model=TauPyModel(model="prem")

file_names=[]
labels=[]
n_list=[872,1482,2405,3326,3941,4534,4922,5405,6026,8112,11289,12181,13332,15821,15988]
for n in n_list:
    labels.append(str(n))
    file_names.append("_"+str(n))
file_names.append("_glad")
file_names.append("_real")
#file_names.append("_mod")
file_names.append("_mod_3")


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
"""
handles=[]
for i,df in enumerate(file_names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker=".", label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends")
plt.close()
"""

phase_list=['SS_S','SS_Sdiff']
ks=[0,1,2,4]
o_l=[8,10,12,34,16]
file_names2=["_real","_8112","_mod_3"]
labels=["Observed","S40RTS","mPREM"]


df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]

df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<140)]
df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
for i,df in enumerate(file_names):
    for w in [0,1]:
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]


file_names=[]
labels=[]
n_list=[872,1482,2405,3326,3941,4534,4922,5405,6026,8112,11289,12181,13332,15821,15988]
for n in n_list:
    labels.append(str(n))
    file_names.append("_"+str(n))
file_names.append("_glad")
file_names.append("_real")
file_names.append("_mod_3") 

for k in range(0,1):
    k=ks[k]
    df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
    df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
   #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))] 
    if k==1:
       df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<110)]
    elif k==0:
         df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]
    elif k==2:
         df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]#df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]

   
    df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    (df_plot["1"].str.split(".").str[0].isin(gsn))]
    df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
    print(len(df_plot))
    
    # labels=["PREM","mPREM","S40RTS"]
    # file_names2=["_prem","_mod_3","_8112"]
    # for i,df in enumerate(file_names2):
    #     plt.figure(1,figsize=(15,10))
    #     # if k==0:
    #     plt.scatter(df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df],np.log(df_plot["8"+df]),alpha=0.5,marker="*",color="black")
    #     # else:
    #     #     plt.scatter(df_plot["7"],df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df],alpha=0.5,marker="*",color="black")
        

    #     plt.grid()
    #     plt.ylabel("SS-S Amplitude Ratios")
    #     plt.xlabel("SS-S CC Differnetial Traveltime")
    #     plt.savefig(plot_dir+"/relative/scatters/cc_amp_"+file_names2[i]+"_"+str(k)+"_scatter.png",dpi=600)
    #     plt.close()

    
    pool = mp.Pool(7)
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

    
    
    
    file_names2=["_real","_8112","_glad"]
    labels=["mPREM","PREM","GLAD-M25"]


    df_plot[str(14+1*5+4)+"_mp_s40"]=df_plot[str(14+1*5+4)+"_8112"]-df_plot[str(14+1*5+4)+"_mprem"]
    df_plot[str(14+0*5+4)+"_mp_s40"]=df_plot[str(14+0*5+4)+"_8112"]-df_plot[str(14+0*5+4)+"_mprem"]

    df_plot[str(14+1*5+4)+"_mp_glad"]=df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+1*5+4)+"_mprem"]
    df_plot[str(14+0*5+4)+"_mp_glad"]=df_plot[str(14+0*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_mprem"]
    fig, [ax0, ax] = plt.subplots(1, 2, figsize=(12.5, 8.6), gridspec_kw={'width_ratios': [1, 2.5]})
    for i,df in enumerate(["_mod_3","_real"]):

        
        y,edges=np.histogram(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df],bins=15,range=(-2,2))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            ax.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3)
            p=np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df])
        elif i==1:
            ax.stairs(y/np.sum(y), edges, color="black",label=labels[i],ls="--",fill=False,lw=3)
            mps=np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df])
        else:
            ax.stairs(y/np.sum(y), edges, color="blue",label=labels[i],ls="--",fill=False,lw=3)
            mps40=np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df])
    ax.axvline(0,color="grey",lw=1)
    #ax.axvline(mps,color="black",ls="--",lw=1) 
    #ax.axvline(mps40,color="blue",ls="--",lw=1)  
    if k==0:
       ax.set_xlabel(r"$\delta \rm lnA^{ss/s}$")
    else:
        ax.set_xlabel(r"$\delta \rm lnA^{ss/s}$")
    ax.set_ylabel("Normalized Frequency")
    
   
    
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.grid(True, which='both', linestyle='dotted', linewidth=1)
    ax.set_xlim(-2,2)
    # fig.tight_layout()
    # plt.savefig(plot_dir+"/relative/histograms/dlnA_mp_"+str(k)+".pdf",dpi=600)
    # plt.close() 

p_depths=np.arange(2500000,2891000,3000)
x=np.zeros(len(p_depths))
prem_vp=np.zeros(len(p_depths))
prem_vs=np.zeros(len(p_depths))
for i,k in enumerate(x):
    x[i]=6371000-p_depths[i]
    if x[i]>3480000 and x[i]<=3630000:
       prem_vp[i] = 15.3891 - 5.3181*(x[i]/6371000.0) + 5.5242*(x[i]/6371000.0)**2 - 2.5514*(x[i]/6371000.0)**3
       prem_vs[i] = 6.9254 + 1.4672*(x[i]/6371000.0) - 2.0834*(x[i]/6371000.0)**2 + 0.9783*(x[i]/6371000.0)**3
    else:
         prem_vp[i] = 24.9520- 40.4673*(x[i]/6371000.0) + 51.4832*(x[i]/6371000.0)**2 - 26.6419*(x[i]/6371000.0)**3
         prem_vs[i] = 11.1671 - 13.7818*(x[i]/6371000.0) + 17.4575*(x[i]/6371000.0)**2 - 9.2777*(x[i]/6371000.0)**3

m_prem_vs=np.zeros(len(p_depths))
for i,k in enumerate(prem_vs):
    if x[i]>3480000 and x[i]<=3671000:
       m_prem_vs[i]=prem_vs[i]*(1-(1-(1/191000)*(x[i]-3480000))/100)
    else:
        m_prem_vs[i]=prem_vs[i]
       #print("done")
           
ax0.text(-0.45,1,"a)",horizontalalignment='center',verticalalignment='center',transform=ax0.transAxes,fontsize=28)
ax.text(-0.1,1,"b)",horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=28)    

ax0.plot(prem_vs,p_depths/1000,label="PREM",color="black",lw=3,linestyle="--")
ax0.plot(m_prem_vs,p_depths/1000,label="mPREM",color="black",lw=3)
ax0.set_ylim(2891,2500)
ax.legend(loc="upper left", fontsize=24)
#ax0.legend(loc="upper right", fontsize=24)
ax0.grid()
fig.tight_layout()
#ax0.legend(loc="upper right")
ax0.set_ylabel("Depth(km)")
ax0.set_xlabel(r"Vs(km/s)")
plt.savefig("panel_new.pdf",dpi=600)
plt.close()
    # for i,df in enumerate(["_mod_3","_mp_s40","_mp_glad"]):
    #     fig=pygmt.Figure()
    #     if k==2:
    #         pygmt.makecpt(cmap='/home/ayon/WORK/plots_final/col_am.cpt',series=[-1.8,1.8],reverse=True,continuous=False) 
    #     else:
    #         pygmt.makecpt(cmap='/home/ayon/WORK/plots_final/col_am.cpt',series=[-1.8,1.8],reverse=True,continuous=False)
    #     df_plot["plot"] = df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df]
    #     df_plot["plot"] = df_plot["plot"].apply(lambda n: -1.8 if n > 1.8  else -1.8 if n < -1.8 else n)
    #     fig.basemap(region="d",projection="R-150/12c",frame="f")
    #     fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.4p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
    #     fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.4p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
    #     fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.25c",fill=df_plot["plot"],cmap=True, transparency=40)
    #     fig.coast(shorelines="1/0.75p",frame="f")
    #     #if i==0:
    #     #   fig.colorbar()
    #     fig.savefig(plot_dir+"/relative/maps/amp_2_"+str(k)+"_"+str(df)+"_receiver_1.5.pdf",dpi=600)
    #     plt.close()

    