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
with open("./../event_list","r") as txt:
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



us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"]         

file_names      = ["_real","_1D","_3D","_3D_2","_3D_15","_prem_3D","_glad","_synt",]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['crimson','indigo','dodgerblue','teal','palegreen','peru','indianred','grey','black']    
lss   =[1,1,1,1,1,1,1,1]
style =['None','solid',"solid","solid","solid","solid",'solid']
labels=["real_data",'SPREMc','S8000',"S11000","S5000","PREMc","GLAD_M25"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None"] 
markers=[".","+","o","*","x","^","v","8",">","<","6"]
alp   =[0.5,1,1,1,1,1,1]
df=[]
plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/figures"
model=TauPyModel(model="prem")


labels=["PREM_3DC","S40RTS","S5000_3DC","S11000_3DC","S40RTS_1DC","GLAD_M25","PREM","PREM_QL6","PREM_QRFSI12","_1D_REF","S40RTS_3DC_3DQ"]
file_names=["_prem_3D","_3D","_3D_5","_3D_11","_1D","_glad","_synt","_prem_16","_3D_atten","_1D_ref","_S80_3D"]

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

phase_list=['S','SS','SSS','ScS','Sdiff']
for k,ph in enumerate(phase_list):
    
    #df_plot=df_plot[moho]
    #df_plot=df_plot[(np.abs(df_plot[str(16+k*7+2)+"_real"])<15) & (np.abs(df_plot[str(16+k*7+2)+"_1D"])<15) &(np.abs(df_plot[str(16+k*7+2)+"_3D"])<15) & (np.abs(df_plot[str(16+k*7+2)+"_prem_3D"])<15)]
    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all_real.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))]
    df_plot=df_plot[(df_plot["7"]>30)]
    #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(file_names):
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+df])<20)]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+1)+df])>0.7)]
        df_plot=df_plot[df_plot[str(18+k*11+5)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+4)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+3)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+3)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+4)+df].notna()]
        df_plot=df_plot[(np.abs(np.log(df_plot[str(18+k*11+5)+df])<1.5))]
        df_plot=df_plot[df_plot[str(18+k*11+6)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+6)+df])<1.5)]
        
    df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    #df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]
    print(len(df_plot["1"]))
    df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    (df_plot["1"].str.split(".").str[0].isin(gsn))]
    print(len(df_plot["1"]))

    
    if k==0:    
        df_plot = df_plot[((df_plot["7"]>30) & (df_plot["7"]<75)) | ((df_plot["7"]>85) & (df_plot["7"]<100))]
        ranges  = [30,60,100]
       
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
        
    
    elif k==1:
        df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
        ranges  = [50,80,110,140]

        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
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
    elif k==2:
        df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>115) & (df_plot["7"]<150))]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        ranges  = [70,100,125,150]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
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
    elif k==3:
        df_plot = df_plot[((df_plot["7"]>5) & (df_plot["7"]<25)) | ((df_plot["7"]>50) & (df_plot["7"]<65)) ]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        ranges=[5,65]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
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
    elif k==4:
        df_plot = df_plot[(df_plot["7"]>100) & (df_plot["7"]<150)]
        ranges  = [100,125,150]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
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

    
    
    
    for i,df in enumerate(file_names):
        fig=plt.figure(figsize=(7.08,3.54))
        ax=fig.add_subplot(1,1,1)
        
        ax.scatter(df_plot["7"],df_plot[str(18+k*11+6)+df],label=labels[i],marker=".",color='blue',s=3)
        

        # Estimate the density of the data
        
        ax.set_ylim(-1.4,1.4)
        ax.text(max(df_plot["7"])-5,0.8,ph)
        ax.grid()
        #plt.tight_layout()
        ax.set_xlabel("Epicentral Distances (Degrees)")
        ax.set_ylabel("$\\chi_{amp}$")
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.95)
        plt.savefig(plot_dir+"/scatter/amp_misfit_"+file_names[i]+"_"+ph+".png")
        plt.close()

        fig=plt.figure(figsize=(7.08,3.54))
        ax=fig.add_subplot(1,1,1)
        
        ax.scatter(df_plot["7"],df_plot[str(18+k*11+3)+df],label=labels[i],marker=".",color='blue',s=3)
        # Estimate the density of the data
        
        ax.set_ylim(-1.4,1.4)
        ax.text(max(df_plot["7"])-5,0.8,ph)
        ax.grid()
        #plt.tight_layout()
        ax.set_xlabel("Epicentral Distances (Degrees)")
        ax.set_ylabel("$\\chi_{env}$")
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.95)
        plt.savefig(plot_dir+"/scatter/env_misfit_"+file_names[i]+"_"+ph+".png")
        plt.close()

        fig=plt.figure(figsize=(7.08,3.54))
        ax=fig.add_subplot(1,1,1)
       

        # Add a density plot
        
        ax.scatter(df_plot["7"],np.log(df_plot[str(18+k*11+5)+df]),label=labels[i],marker=".",color='blue',s=3)
        # Estimate the density of the data
        
        ax.set_ylim(-1.4,1.4)
        ax.text(max(df_plot["7"])-5,0.8,ph)
        ax.grid()
        #plt.tight_layout()
        ax.set_xlabel("Epicentral Distances (Degrees)")
        ax.set_ylabel("$\\chi_{amp_2}$")
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.95)
        plt.savefig(plot_dir+"/scatter/amp_2_"+file_names[i]+"_"+ph+".png")
        plt.close()
 
        #print(df_plot)
        for q in range(len(ranges)):
            if q==0:
                df_plots=df_plot
            else:
                 df_plots = df_plot[((df_plot["7"]>ranges[q-1]) & (df_plot["7"]<ranges[q])) ]
        ### Plotting Global Maps for Amplitude Ratio
            fig=pygmt.Figure()
            if k==2:
               pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4]) 
            else:
                pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4])
            fig.basemap(region="d",projection="N-150/12c",frame=True)
            
            
            fig.plot(data=df_plots[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
            fig.plot(data=df_plots[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
            fig.plot(x=df_plots.mid_lon,y=df_plots.mid_lat,style="c0.2c",fill=np.log(df_plots[str(18+k*11+5)+df]),pen="0.000001p",cmap=True,transparency=40)
            fig.coast(shorelines=True, frame=True)
            #if i==0:
            #   fig.colorbar()
            fig.savefig(plot_dir+"/global/amp_"+str(q)+"_"+file_names[i]+"_"+ph+".png",dpi=600)
            plt.close()

        ### Plotting Global Maps for Envelope Ratio
            fig=pygmt.Figure()

            if k==2:
               pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4]) 
            else:
                pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4])
            fig.basemap(region="d",projection="N-150/12c",frame=True)
            #fig.coast(shorelines=True, frame=True)   
            fig.plot(data=df_plots[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
            fig.plot(data=df_plots[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)
            #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
            #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
            fig.plot(x=df_plots.mid_lon,y=df_plots.mid_lat,style="c0.2c",fill=df_plots[str(18+k*11+3)+df],pen="0.000001p",cmap=True,transparency=40)
            fig.coast(shorelines=True, frame=True)
            #if i==0:
            #   fig.colorbar()
            fig.savefig(plot_dir+"/global/env_"+str(q)+"_"+file_names[i]+"_"+ph+".png",dpi=600)
            plt.close()
        
            fig=pygmt.Figure()
            if k==2:
               pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4]) 
            else:
                pygmt.makecpt(cmap='polar',reverse=True,series=[-1.4,1.4])
            fig.basemap(region="d",projection="N-150/12c",frame=True)
            #fig.coast(shorelines=True, frame=True)   
            fig.plot(data=df_plots[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
            fig.plot(data=df_plots[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)
            #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
            #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
            fig.plot(x=df_plots.mid_lon,y=df_plots.mid_lat,style="c0.5c",fill=df_plots[str(18+k*11+6)+df],pen="0.000001p",cmap=True,transparency=40)
            fig.coast(shorelines=True, frame=True)
            #if i==0:
            #   fig.colorbar()
            fig.savefig(plot_dir+"/global/amp_2_"+str(q)+"_"+file_names[i]+"_"+ph+".png",dpi=600)
            plt.close()

            
    
    

