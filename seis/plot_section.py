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


ev_list=[]
with open("events","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter=' ')
     for row in reader:
        ev_list.append(row[0])
#ev_list=["090297B"]

def plot_sec(ev):
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

    colors=["blue","green","red","purple","brown","pink","gray","olive","cyan","orange","palegreen"]
    model=TauPyModel(model="prem")
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+ev+".T017-040s.proc_synt.h5",mode="r")
    evla=float(ds.events[0].origins[0].latitude)
    evlo=float(ds.events[0].origins[0].longitude)
    evdp=float(ds.events[0].origins[0].depth)/1000.0
    stations=ds.waveforms.list()
    dist_list=[]
    arr_list=[]
    for u,comp in enumerate(["T","R","Z"]):
        plt.figure(u,figsize=(20,30))
        
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
                        data=ds.waveforms[st].proc_synt.select(component=comp)[0]
                        sampling_rate = data.stats.sampling_rate
                        #print(arrival)
                        arr_ind       = int(arrival*sampling_rate)
                        p_data        = data.data[arr_ind-int(50*sampling_rate):arr_ind+int(1000*sampling_rate)]
                        if len(p_data)!=0:
                           p_data        = p_data/max(np.abs(p_data))
                           time          = np.linspace(-50,1000,len(p_data))
                           plt.plot(time,5*p_data+dist,linewidth=0.5,color="black",alpha=0.7)
                        if u==0:
                           dist_list.append(dist)
                           arr_list.append(arrival)
                    except IndexError:
                        continue
        
        for p,ph in enumerate(["S","sS","SKS","SKKS","Sdiff","S410S","S660S","S^660S","S^440S","SS","sSdiff"]):
            d=[]
            t=[]
            for l,dists in enumerate(dist_list):
                arr=model.get_travel_times(evdp,dists,phase_list=[ph])
                if len(arr)!=0:
                    d.append(dists)
                    t.append(arr[0].time-arr_list[l]-48)
            plt.scatter(t,d,label=ph,color=colors[p],marker="o",s=300,alpha=0.7)

            
                    
        


        plt.legend(loc="lower right")
        plt.xlim(-50,1000)
        plt.xlabel("Time Relative to S/Sdiff arrival")
        plt.ylabel("Epicentral Distance(deg)")
        plt.title("Long="+str(evlo)+" Lat="+str(evla)+" Depth="+str(evdp) +" for " + comp +" component")
        plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/sections/"+ev+"_"+comp+"_"+str(evdp)+"_section.png")
        
        plt.close() 
    return 0 

pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
pool.map(plot_sec, ev_list)
pool.close()
pool.join()          

           


    
