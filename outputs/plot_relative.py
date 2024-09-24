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
#from mpi4py import MPI

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


    
         

us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"]
file_names      = ["_real","_prem_3D","_3D","_3D_2","_3D_15","_1D","_glad","_1D_ref","_prem_q16","_3D_atten","_1D"]
colors=['peru','black','forestgreen','magenta','mediumpurple','deepskyblue','darkorange','grey','crimson']   
colors=["black","darkorange","darkmagenta","royalblue","forestgreen","red"]  
lss   =[2,1.5,1,0.5,1,1,1,1,1,1,1,1,1,1]
style =['solid','solid',"solid","solid","solid","solid","solid","solid","solid"]
labels=["Real Data","PREMc","S8000","S11000","S5000","SPREMc","GLAD_M25","1D_REF","PREM_QL6","PREM_QRFSI12","S40RTS_PREM_Crust"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None","None","None","None","None"] 
markers=["v",".","+","o","*","x","^"]
alp   =[1,1,1,1,1,1,1]
df=[]
plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/figures/relative/"
model=TauPyModel(model="prem")





# labels=["S40RTS_3DC","S40RTS_3DC_QRFSI12","PREM_QL6","Real Data"]
# file_names=["_3D","_S80_3D","_prem_16","_real"]

# labels=["Real Data","GLAD_M25","S40RTS_3DC","S40RTS_1DC","PREM_3DC"]
# file_names=["_real","_glad","_3D","_1D","_prem_3D"]
labels=["PREM_3DC","S40RTS_3DC","S40RTS_1DC","GLAD_M25","Real Data"]
file_names=["_prem_3D","_3D","_1D","_glad","_real"]

all_f=["_real","_glad","_3D","_1D","_prem_3D","_3D_11","_3D_5"]
#for i,f in enumerate(file_names):
#    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt",header=None,delimiter=" ",skipinitialspace=True)
 #   df.append(temp)
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

handles=[]
for i,df in enumerate(file_names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker=markers[i], label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends_2")
plt.close()

"""
SS/S Ratios
"""
df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SS_S_all.txt",skipinitialspace=True,delimiter=",")
#df_plot=df[df[str(8)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]


df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
ranges=[50,75,100]
#df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
for i,df in enumerate(all_f):
    df_plot=df_plot[(np.abs(df_plot[str(18+0*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+0*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(18+1*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+1*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(np.log((df_plot[str(8)+df])))<1.5)]

    
#df_plot=df_plot[(np.abs(df_plot[str(18)+"_real"])<16) & (np.abs(df_plot[str(18)+"_1D"])<16) &(np.abs(df_plot[str(18)+"_3D"])<16) & (np.abs(df_plot[str(18)+"_prem_3D"])<16)]
#df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<16) & (np.abs(df_plot[str(24)+"_1D"])<16) &(np.abs(df_plot[str(24)+"_3D"])<16) & (np.abs(df_plot[str(24)+"_prem_3D"])<16)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
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

ls=["7"]
df_plot_2=df_plot.copy()


plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    ls.append(str(8)+df)
print(ls)
df_dist = df_plot_2[ls].copy()
for cl in ls[1:]:
    df_dist[cl] = np.log(df_plot_2[cl])
df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
df_dist_2["dist"] = np.arange(32.5,150,2.5)
print(df_dist)
for i,cl in enumerate(ls[1:]):
    plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5,markersize=4)
plt.ylim(-0.5,0.5)
plt.xlabel("Epicentral Distance")
plt.ylabel("SS/S ratio")
plt.tight_layout()
plt.savefig(plot_dir+"/histograms/ep_SS_S.png",dpi=600)
plt.close()

for i,df in enumerate(file_names):

        
    for q in range(len(ranges)):
        if q==0:
                df_plots=df_plot
        else:
                df_plots = df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ] 
        ### Plotting Global Maps for Amplitude Ratio
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        
        
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.2c",fill=np.log(df_plot[str(8)+df]),pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        #if i==0:
        #    fig.colorbar()
        fig.savefig(plot_dir+"/global/amp_"+file_names[i]+"_"+str(q)+"_SS_S.png",dpi=600)
        plt.close()

for q in range(len(ranges)):
    means=[]
    std=[]
    if q==0:
        df_plots=df_plot
    else:
        df_plots=df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ]
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
        y, edges = np.histogram(np.log(df_plot[str(8)+df]), bins=14, range=(-1.4,1.4))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)
        means.append(np.mean(np.log(df_plot[str(8)+df])))
        std.append(np.std(np.log(df_plot[str(8)+df])))
    plt.text(1.2,0.8*max(y/np.sum(y)),"SS/S")
    plt.xlim(-1.4,1.4)
    #plt.legend(fontsize=12)
    plt.xlabel("$\\chi_{SS/S}$")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/histograms/line_SS_S_"+str(q)+".png",dpi=600)
    plt.close()
    with open(plot_dir+"/outputs_relative.txt","a") as txt:
        txt.write("SS/S Ratio mean ")
        for m,mea in enumerate(means):
            txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" ")
        txt.write(str(np.sum(y))+" "+"\n")

for i,df in enumerate(file_names[0:4]):
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
            
    ax.plot(np.linspace(-1.4,1.4,100), np.linspace(-1.4,1.4,100), linestyle='--', color='black')
    ax.plot(np.linspace(-1.4,1.4,100), -np.linspace(-1.4,1.4,100), linestyle='--', color='black')

    
    ax.scatter(np.log(df_plot[str(8)+"_real"]),np.log(df_plot[str(8)+df]),marker=".",alpha=0.6)
    ax.set_xlabel("$\\chi_{SS/S}$ for Real Data")
    ax.set_ylabel("$\\chi_{SS/S}$ for "+labels[i])
    cc=df_plot[str(8)+"_real"].corr(df_plot[str(8)+df])
    ax.set_title("cc = "+str(cc))
    ax.set_xlim([-1.4, 1.4])
    ax.set_ylim([-1.4, 1.4])
    plt.tight_layout()
    plt.savefig(plot_dir+"/scatter/Scatter_SS_S"+df+".png",dpi=600)
    plt.close()
misfits_SS_S=[]
i_s=[]
for i,df in enumerate(file_names[0:4]):
    print(df)
    i_s.append(i)
    mfs=(np.log(df_plot[str(8)+"_real"]/df_plot[str(8)+df])**2).mean()
    misfits_SS_S.append(mfs)


"""
SS/Sdiff Ratios
"""
df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SS_Sdiff_all.txt",skipinitialspace=True,delimiter=",")
#df_plot=df[df[str(8)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]

# df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
#     ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
#     ~(df_plot["1"].str.split(".").str[0].isin(us)) &
#     ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
#     (df_plot["1"].str.split(".").str[0].isin(gsn))]

df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
ranges=[100,125,140]
#df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
for i,df in enumerate(all_f):
    df_plot=df_plot[(np.abs(df_plot[str(18+4*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+4*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(18+1*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+1*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(np.log((df_plot[str(10)+df])))<1.5)]

    
#df_plot=df_plot[(np.abs(df_plot[str(18)+"_real"])<15) & (np.abs(df_plot[str(18)+"_1D"])<15) &(np.abs(df_plot[str(18)+"_3D"])<15) & (np.abs(df_plot[str(18)+"_prem_3D"])<15)]
#df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<15) & (np.abs(df_plot[str(24)+"_1D"])<15) &(np.abs(df_plot[str(24)+"_3D"])<15) & (np.abs(df_plot[str(24)+"_prem_3D"])<15)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
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

ls=["7"]
df_plot_2=df_plot.copy()


plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    ls.append(str(10)+df)
print(ls)
df_dist = df_plot_2[ls].copy()
for cl in ls[1:]:
    df_dist[cl] = np.log(df_plot_2[cl])
df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
df_dist_2["dist"] = np.arange(32.5,150,2.5)
print(df_dist)
for i,cl in enumerate(ls[1:]):
    plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5,markersize=4)
plt.ylim(-0.5,0.5)
plt.xlabel("Epicentral Distance")
plt.ylabel("SS/Sdiff ratio")
plt.tight_layout()
plt.savefig(plot_dir+"/histograms/ep_SS_Sdiff.png",dpi=600)
plt.close()

for i,df in enumerate(file_names):

    for q in range(len(ranges)):
        if q==0:
                df_plots=df_plot
        else:
                df_plots = df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ]   

        ### Plotting Global Maps for Amplitude Ratio
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        
        
        fig.plot(data=df_plots[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plots[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
        fig.plot(x=df_plots.mid_lon,y=df_plots.mid_lat,style="c0.2c",fill=np.log(df_plots[str(10)+df]),pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        #if i==0:
        #    fig.colorbar()
        fig.savefig(plot_dir+"/global/amp_"+file_names[i]+"_"+str(q)+"_SS_Sdiff.png",dpi=600)
        plt.close()


for q in range(len(ranges)):
    means=[]
    std=[]
    if q==0:
        df_plots=df_plot
    else:
        df_plots=df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ]
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
        y, edges = np.histogram(np.log(df_plots[str(10)+df]), bins=14, range=(-1.4,1.4))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)
        means.append(np.mean(np.log(df_plots[str(10)+df])))
        std.append(np.std(np.log(df_plots[str(10)+df])))
    plt.text(1.2,0.8*max(y/np.sum(y)),"SS/Sdiff")
    plt.xlim(-1.4,1.4)
    #plt.legend(fontsize=12)
    plt.xlabel("$\\chi_{SS/Sdiff}$")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/histograms/line_SS_Sdiff_"+str(q)+".png",dpi=600)
    plt.close()

    with open(plot_dir+"/outputs_relative.txt","a") as txt:
            txt.write("SS/Sdiff Ratio mean "+"\n")
            for m,mea in enumerate(means):
                txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" for range "+str(q)+"\n")
            txt.write(str(np.sum(y))+" "+"\n")
for i,df in enumerate(file_names[0:4]):
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
            
    ax.plot(np.linspace(-1.4,1.4,100), np.linspace(-1.4,1.4,100), linestyle='--', color='black')
    ax.plot(np.linspace(-1.4,1.4,100), -np.linspace(-1.4,1.4,100), linestyle='--', color='black')

    
    ax.scatter(np.log(df_plot[str(10)+"_real"]),np.log(df_plot[str(10)+df]),marker=".",alpha=0.6)
    ax.set_xlabel("$\\chi_{SS/Sdiff}$ for Real Data")
    ax.set_ylabel("$\\chi_{SS/Sdiff}$ for "+labels[i])
    cc=df_plot[str(10)+"_real"].corr(df_plot[str(10)+df])
    ax.set_title("cc = "+str(cc))
    ax.set_xlim([-1.4, 1.4])
    ax.set_ylim([-1.4, 1.4])
    plt.tight_layout()
    plt.savefig(plot_dir+"/scatter/Scatter_SS_Sdiff"+df+".png",dpi=600)
    plt.close()
misfits_SS_Sdiff=[]
i_s=[]
for i,df in enumerate(file_names[0:4]):
    
    i_s.append(i)
    mfs=(np.log(df_plot[str(10)+"_real"]/df_plot[str(10)+df])**2).mean()
    misfits_SS_Sdiff.append(mfs)

"""
SSS/SS Ratios 
"""
df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS_SS_all.txt",skipinitialspace=True,delimiter=",")
#df_plot=df[df[str(10)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
# df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
#     ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
#     ~(df_plot["1"].str.split(".").str[0].isin(us)) &
#     ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
#     (df_plot["1"].str.split(".").str[0].isin(gsn))]

df_plot = df_plot[(((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<140)))]
ranges=[70,100,125,140]
#df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
for i,df in enumerate(all_f):
    df_plot=df_plot[(np.abs(df_plot[str(18+2*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+2*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(18+1*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+1*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(np.log((df_plot[str(12)+df])))<1.5)]

    
#df_plot=df_plot[(np.abs(df_plot[str(18)+"_real"])<15) & (np.abs(df_plot[str(18)+"_1D"])<15) &(np.abs(df_plot[str(18)+"_3D"])<15) & (np.abs(df_plot[str(18)+"_prem_3D"])<15)]
#df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<15) & (np.abs(df_plot[str(24)+"_1D"])<15) &(np.abs(df_plot[str(24)+"_3D"])<15) & (np.abs(df_plot[str(24)+"_prem_3D"])<15)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
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

ls=["7"]
df_plot_2=df_plot.copy()


plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    ls.append(str(12)+df)
print(ls)
df_dist = df_plot_2[ls].copy()
for cl in ls[1:]:
    df_dist[cl] = np.log(df_plot_2[cl])
df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
df_dist_2["dist"] = np.arange(32.5,150,2.5)
print(df_dist)
for i,cl in enumerate(ls[1:]):
    plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5,markersize=4)
plt.ylim(-0.5,0.5)
plt.xlabel("Epicentral Distance")
plt.ylabel("SSS/SS ratio")
plt.tight_layout()
plt.savefig(plot_dir+"/histograms/ep_SSS_SS.png",dpi=600)
plt.close()

for i,df in enumerate(file_names):

        

        ### Plotting Global Maps for Amplitude Ratio
    for q in range(len(ranges)):
        if q==0:
                df_plots=df_plot
        else:
                df_plots = df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ] 
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        
        
        fig.plot(data=df_plots[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plots[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
        fig.plot(x=df_plots.mid_lon,y=df_plots.mid_lat,style="c0.2c",fill=np.log(df_plots[str(12)+df]),pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        #if i==0:
        #    fig.colorbar()
        fig.savefig(plot_dir+"/global/global_amp_"+file_names[i]+"_"+str(q)+"_SSS_SS.png",dpi=600)
        plt.close()
for q in range(len(ranges)):
    means=[]
    std=[]
    if q==0:
        df_plots=df_plot
    else:
        df_plots=df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ]
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
        y, edges = np.histogram(np.log(df_plots[str(12)+df]), bins=14, range=(-1.4,1.4))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)
        means.append(np.mean(np.log(df_plots[str(12)+df])))
        std.append(np.std(np.log(df_plots[str(12)+df])))
    plt.text(1.2,0.8*max(y/np.sum(y)),"SSS/SS")
    plt.xlim(-1.4,1.4)
    #plt.legend(fontsize=12)
    plt.xlabel("$\\chi_{SSS/SS}$")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/histograms/line_SSS_SS"+"_"+str(q)+".png",dpi=600)
    plt.close()

    with open(plot_dir+"/outputs_relative.txt","a") as txt:
            txt.write("SSS/SS Ratio mean "+"\n")
            for m,mea in enumerate(means):
                txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" for range "+str(q)+"\n")
            txt.write(str(np.sum(y))+" "+"\n")

for i,df in enumerate(file_names[0:4]):
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
            
    ax.plot(np.linspace(-1.4,1.4,100), np.linspace(-1.4,1.4,100), linestyle='--', color='black')
    ax.plot(np.linspace(-1.4,1.4,100), -np.linspace(-1.4,1.4,100), linestyle='--', color='black')

    
    ax.scatter(np.log(df_plot[str(12)+"_real"]),np.log(df_plot[str(12)+df]),marker=".",alpha=0.6)
    ax.set_xlabel("$\\chi_{SSS/SS}$ for Real Data")
    ax.set_ylabel("$\\chi_{SSS/SS}$ for "+labels[i])
    cc=df_plot[str(12)+"_real"].corr(df_plot[str(12)+df])
    ax.set_title("cc = "+str(cc))
    ax.set_xlim([-1.4, 1.4])
    ax.set_ylim([-1.4, 1.4])
    plt.tight_layout()
    plt.savefig(plot_dir+"/scatter/Scatter_SSS_SS"+df+".png",dpi=600)
    plt.close()
misfits_SSS_SS=[]
i_s=[]
for i,df in enumerate(file_names[0:4]):
    
    i_s.append(i)
    mfs=(np.log(df_plot[str(12)+"_real"]/df_plot[str(12)+df])**2).mean()
    misfits_SSS_SS.append(mfs)

"""
SSS/Sdiff Ratios 
"""
df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS_Sdiff_all.txt",skipinitialspace=True,delimiter=",")
#df_plot=df[df[str(10)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
# df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
#     ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
#     ~(df_plot["1"].str.split(".").str[0].isin(us)) &
#     ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
#     (df_plot["1"].str.split(".").str[0].isin(gsn))]

df_plot = df_plot[(((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<150)))]
ranges=[100,125,150]
#df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
for i,df in enumerate(all_f):
    df_plot=df_plot[(np.abs(df_plot[str(18+2*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+2*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(df_plot[str(18+4*11+2)+df])<20)]
    df_plot=df_plot[(np.abs(df_plot[str(18+4*11+1)+df])>0.7)]
    df_plot=df_plot[(np.abs(np.log((df_plot[str(16)+df])))<1.5)]

    
#df_plot=df_plot[(np.abs(df_plot[str(18)+"_real"])<15) & (np.abs(df_plot[str(18)+"_1D"])<15) &(np.abs(df_plot[str(18)+"_3D"])<15) & (np.abs(df_plot[str(18)+"_prem_3D"])<15)]
#df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<15) & (np.abs(df_plot[str(24)+"_1D"])<15) &(np.abs(df_plot[str(24)+"_3D"])<15) & (np.abs(df_plot[str(24)+"_prem_3D"])<15)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
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

ls=["7"]
df_plot_2=df_plot.copy()


plt.figure(1,figsize=(7.08,3.54))
for i,df in enumerate(file_names):
    ls.append(str(16)+df)
print(ls)
df_dist = df_plot_2[ls].copy()
for cl in ls[1:]:
    df_dist[cl] = np.log(df_plot_2[cl])
df_dist_2 = df_dist.groupby(pd.cut(df_dist["7"],np.arange(30,150,2.5))).mean()
df_dist_2["dist"] = np.arange(32.5,150,2.5)
print(df_dist)
for i,cl in enumerate(ls[1:]):
    plt.plot(df_dist_2["dist"], df_dist_2[cl],label=labels[i],marker=markers[i],color=colors[i],lw=0.5,markersize=4)
plt.ylim(-0.5,0.5)
plt.xlabel("Epicentral Distance")
plt.ylabel("SSS/Sdiff ratio")
plt.tight_layout()
plt.savefig(plot_dir+"/histograms/ep_SSS_Sdiff.png",dpi=600)
plt.close()

for i,df in enumerate(file_names):

        

        ### Plotting Global Maps for Amplitude Ratio
    for q in range(len(ranges)):
        if q==0:
                df_plots=df_plot
        else:
                df_plots = df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ] 
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        
        
        fig.plot(data=df_plots[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plots[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
        fig.plot(x=df_plots.mid_lon,y=df_plots.mid_lat,style="c0.2c",fill=np.log(df_plots[str(16)+df]),pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        #if i==0:
        #    fig.colorbar()
        fig.savefig(plot_dir+"/global/global_amp_"+file_names[i]+"_"+str(q)+"_SSS_Sdiff.png",dpi=600)
        plt.close()
for q in range(len(ranges)):
    means=[]
    std=[]
    if q==0:
        df_plots=df_plot
    else:
        df_plots=df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ]

    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
        y, edges = np.histogram(np.log(df_plots[str(16)+df]), bins=14, range=(-1.4,1.4))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i],markersize=4)
        means.append(np.mean(np.log(df_plots[str(16)+df])))
        std.append(np.std(np.log(df_plots[str(16)+df])))
    plt.text(1.2,0.8*max(y/np.sum(y)),"SSS/Sdiff")
    plt.xlim(-1.4,1.4)
    #plt.legend(fontsize=12)
    plt.xlabel("$\\chi_{SSS/Sdiff}$")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/histograms/line_SSS_Sdiff"+"_"+str(q)+".png",dpi=600)
    plt.close()

    with open(plot_dir+"/outputs_relative.txt","a") as txt:
            txt.write("SSS/Sdiff Ratio mean "+"\n")
            for m,mea in enumerate(means):
                txt.write(str(np.exp(mea))+" SD "+str(np.exp(std[m]))+" "+str(np.exp(-std[m]))+" for range "+str(q))
            txt.write(str(np.sum(y))+" "+"\n")

for i,df in enumerate(file_names[0:4]):
    fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
            
    ax.plot(np.linspace(-1.4,1.4,100), np.linspace(-1.4,1.4,100), linestyle='--', color='black')
    ax.plot(np.linspace(-1.4,1.4,100), -np.linspace(-1.4,1.4,100), linestyle='--', color='black')

    
    ax.scatter(np.log(df_plot[str(16)+"_real"]),np.log(df_plot[str(16)+df]),marker=".",alpha=0.6)
    ax.set_xlabel("$\\chi_{SSS/Sdiff}$ for Real Data")
    ax.set_ylabel("$\\chi_{SSS/Sdiff}$ for "+labels[i])
    cc=df_plot[str(16)+"_real"].corr(df_plot[str(16)+df])
    ax.set_title("cc = "+str(cc))
    ax.set_xlim([-1.4, 1.4])
    ax.set_ylim([-1.4, 1.4])
    plt.tight_layout()
    plt.savefig(plot_dir+"/scatter/Scatter_SSS_Sdiff"+df+".png",dpi=600)
    plt.close()
misfits_SSS_Sdiff=[]
i_s=[]
for i,df in enumerate(file_names[0:4]):
    
    i_s.append(i)
    mfs=(np.log(df_plot[str(16)+"_real"]/df_plot[str(16)+df])**2).mean()
    misfits_SSS_Sdiff.append(mfs)


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


fig,ax=plt.subplots(figsize=(7.08,4))
ax.plot(i_s,misfits_SS_S,linestyle="--",marker=markers[0],color=colors[0],label="$\\chi_{SS/S}$",lw=1)
ax.plot(i_s,misfits_SS_Sdiff,linestyle="--",marker=markers[1],color=colors[1],label="$\\chi_{SS/Sdiff}$",lw=1)
ax.plot(i_s,misfits_SSS_SS,linestyle="--",marker=markers[2],color=colors[2],label="$\\chi_{SSS/SS}$",lw=1)
ax.plot(i_s,misfits_SSS_Sdiff,linestyle="--",marker=markers[3],color=colors[3],label="$\\chi_{SSS/Sdiff}$",lw=1)
ax.set_xticks(i_s)
ax.set_xticklabels(labels[0:4])
ax.set_xlabel("Global Models")
ax.set_ylabel("Average Misfit")
ax.legend()
ax.grid() 
plt.savefig(plot_dir+"/histograms/misfit_relative.png")