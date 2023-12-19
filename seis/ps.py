import matplotlib.pyplot as plt
import pandas as pd
from pyasdf import ASDFDataSet
from scipy.signal import hilbert
import multiprocessing as mp
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
import csv
from obspy.taup import TauPyModel

SMALL_SIZE = 20
MEDIUM_SIZE = 24
BIGGER_SIZE = 26

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)

ev_list=[]
with open("events","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter=' ')
     for row in reader:
        ev_list.append(row[0])
ev_list=["032899F"]
model=TauPyModel(model="prem")
for ev in ev_list:
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+ev+".T017-040s.proc_synt.h5",mode="r")
    evla=float(ds.events[0].origins[0].latitude)
    evlo=float(ds.events[0].origins[0].longitude)
    evdp=float(ds.events[0].origins[0].depth)/1000.0
    stations=ds.waveforms.list()
    plt.figure(1,figsize=(20,30))
    for st in stations:
        try:
            stla=ds.waveforms[st].coordinates['latitude']
            stlo=ds.waveforms[st].coordinates['longitude']
        except KeyError:
               stla=np.nan
               stlo=np.nan
        if stla==stla and stlo==stlo:
            k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['S'])
            ph=0
            if len(k) != 0:
                dist=k[0].distance
                arrival=k[0].time
                ph=0
            else:
                k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['Sdiff']) 
                if len(k) != 0:
                    dist=k[0].distance
                    arrival=k[0].time
                    ph=1
                else:
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['ScS'])
                    if len(k) != 0:
                        dist=k[0].distance
                        arrival=k[0].time
                        ph=2
                    else:
                        dist=np.nan
                        arrival=np.nan
            if dist is not np.nan and arrival is not np.nan:
                try:
                  data=ds.waveforms[st].proc_synt.select(component="T")[0]
                  sampling_rate = data.stats.sampling_rate
                  print(arrival)
                  arr_ind       = int(arrival*sampling_rate)
                  p_data        = data.data[arr_ind-200:arr_ind+4000]
                  p_data        = p_data/max(p_data)
                  time          = np.linspace(-50,1000,len(p_data))
                  plt.plot(time,5*p_data+dist,linewidth=0.5,color="black")
                except [IndexError,ValueError]:
                       continue
                
                if ph==0:
                    plt.scatter(arrival-arrival-50,dist,marker="|",color="red",s=30,label="S Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['SS'])
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="blue",s=30,label="SS Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['SSS'])
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="gold",s=30,label="SSS Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['ScS']) 
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="skyblue",s=30,label="ScS Phase") 
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['Sdiff']) 
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="green",s=30,label="Sdiff Phase")
                if ph==1:
                    plt.scatter(arrival-arrival-50,dist,marker="|",color="green",s=30,label="Sdiff Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['SS'])
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="blue",s=30,label="SS Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['SSS'])
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="gold",s=30,label="SSS Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['ScS']) 
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="skyblue",s=30,label="ScS Phase") 
                if ph==2:
                    plt.scatter(arrival-arrival-50,dist,marker="|",color="skyblue",s=30,label="ScS Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['SS'])
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="blue",s=30,label="SS Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['SSS'])
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="gold",s=30,label="SSS Phase")
                    k=model.get_ray_paths_geo(evdp,evla,evlo,stla,stlo,['Sdiff']) 
                    if len(k)!=0:
                        p_arr=k[0].time
                        plt.scatter(p_arr-arrival-50,dist,marker="|",color="green",s=30,label="Sdiff Phase") 
                
    



    plt.xlabel("Time Relative to S/Sdiff arrival")
    plt.ylabel("Epicentral Distance(deg)")
    plt.title("Long="+str(evlo)+" Lat="+str(evla)+" Depth="+str(evdp))
    plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/sections/"+ev+"_section.png")
    
    plt.close()            


           


    
