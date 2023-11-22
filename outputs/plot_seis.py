import matplotlib.pyplot as plt
import pandas as pd
from pyasdf import ASDFDataSet
from scipy.signal import hilbert
import multiprocessing as mp
import numpy as np
import os
from mpl_toolkits.basemap import Basemap

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/S_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[(df["7"]>30) & (df["7"]<70)]

def plot(event,station,time,shift,loc,rat):
    names=["proc_obsd_1D_crust","proc_obsd_3D_crust","proc_obsd_glad","proc_real_data"]
    tags =["proc_obsd_1","proc_obsd_3","proc_obsd_25","proc_real_data"]
    labels=["S40RTS+1D Crust","S40RTS+Crust2.0","GLAD_M25","Real Data"]
    colors=["red","skyblue","green","blueviolet"]
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s.proc_synt.h5",mode="r")
    data=ds.waveforms[station].proc_synt.select(component="T")[0]
    sampling_rate = data.stats.sampling_rate
    dt            = data.stats.delta
    anltc         = hilbert(data.data)
    envp          = np.abs(anltc)
    arr_ind       = int(time*sampling_rate)
    minus         = int(50*sampling_rate)
    plus          = int(50*sampling_rate)
    p_data        = data.data[arr_ind-minus:arr_ind+plus]
    p_envp        = envp[arr_ind-minus:arr_ind+plus]
    time_series          = data.times()[arr_ind-minus:arr_ind+plus]
    elat          = loc[1]
    elon          = loc[0]
    slat          = loc[4]
    slon         = loc[3]
    EVT           = np.array([0,loc[2]])
    STA           = np.array([loc[5],0])

    
    
    fig, (a0, a1) = plt.subplots(2, 1, height_ratios=[2, 1],figsize=(25,15))
    a0.plot(time_series,p_data,color="black",lw=1,label="1D PREM")
    a0.plot(time_series,p_envp,ls="--",color="black",lw=1,label="1D PREM envelope")
    del ds
    
    for j,n in enumerate(names):
        ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s."+n+".h5",mode="r")
        data=ds.waveforms[station][tags[j]].select(component="T")[0]
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        anltc         = hilbert(data.data)
        envp          = np.abs(anltc)
        arr_ind       = int((time+shift[j])*sampling_rate)
        minus         = int(50*sampling_rate)
        plus          = int(50*sampling_rate)
        p_data        = data.data[arr_ind-minus:arr_ind+plus]
        p_envp        = envp[arr_ind-minus:arr_ind+plus]
        #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
        del ds
        a0.plot(time_series,p_data,color=colors[j],lw=1,label=labels[j]+"(shift = "+str(shift[j])+")")
        a0.plot(time_series,p_envp,ls="--",color=colors[j],lw=1,label=labels[j]+" envelopes")
    a0.legend(loc="upper right")
    map = Basemap(projection='robin', lon_0=-150, resolution="c",ax=a1)
    map.drawmapboundary(fill_color='#cccccc')
    map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)
    map.scatter(elon, elat, s=80, zorder=10,color="red", marker="o",edgecolor="k")
    map.scatter(slon, slat, s=80, zorder=10,color="lime", marker="o",edgecolor="k")
    if rat>=0:
       map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='red')
    else:
        map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='blue')
    plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/S_2/"+event+"_"+station+".png")
    plt.close(fig)

#for i,rows in df.iterrows():
#    plot(rows[0],rows[1],rows[8],[rows[16],rows[22],rows[28],rows[22]],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]])
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot, [(rows[0],rows[1],rows[8],[rows[16],rows[22],rows[28],rows[10]],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]],np.log(rows[23])/np.log(rows[17]))  for i,rows in df.iterrows()])
pool.close()
pool.join()

df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SS_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[(df["7"]>50) & (df["7"]<140)]

def plot(event,station,time,shift):
    names=["proc_obsd_1D_crust","proc_obsd_3D_crust","proc_obsd_glad","proc_real_data"]
    tags =["proc_obsd_1","proc_obsd_3","proc_obsd_25","proc_real_data"]
    labels=["S40RTS+1D Crust","S40RTS+Crust2.0","GLAD_M25","Real Data"]
    colors=["red","skyblue","green","blueviolet"]
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s.proc_synt.h5")
    data=ds.waveforms[station].proc_synt.select(component="T")[0]
    sampling_rate = data.stats.sampling_rate
    dt            = data.stats.delta
    anltc         = hilbert(data.data)
    envp          = np.abs(anltc)
    arr_ind       = int(time*sampling_rate)
    minus         = int(50*sampling_rate)
    plus          = int(50*sampling_rate)
    p_data        = data.data[arr_ind-minus:arr_ind+plus]
    p_envp        = envp[arr_ind-minus:arr_ind+plus]
    time_series          = data.times()[arr_ind-minus:arr_ind+plus]
    plt.figure(1,figsize=(15,5))
    plt.plot(time_series,p_data,color="black",linewidth=1,label="1D PREM")
    plt.plot(time_series,p_envp,linestyle="--",color="black",linewidth=1,label="1D PREM envelope")
    del ds
    
    for j,n in enumerate(names):
        ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s."+n+".h5")
        data=ds.waveforms[station][tags[j]].select(component="T")[0]
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        anltc         = hilbert(data.data)
        envp          = np.abs(anltc)
        arr_ind       = int((time+shift[j])*sampling_rate)
        minus         = int(50*sampling_rate)
        plus          = int(50*sampling_rate)
        p_data        = data.data[arr_ind-minus:arr_ind+plus]
        p_envp        = envp[arr_ind-minus:arr_ind+plus]
        #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
        del ds
        plt.plot(time_series,p_data,color=colors[j],linewidth=1,label=labels[j]+"(shift = "+str(shift[j])+")")
        plt.plot(time_series,p_envp,linestyle="--",color=colors[j],linewidth=1,label=labels[j]+" envelopes")
    plt.legend(loc=(1.04, 1c))
    plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SS/"+event+"_"+station+".png")
    plt.close()


pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot, [(rows[0],rows[1],rows[8],[rows[9],rows[15],rows[21],rows[27]])  for i,rows in df.iterrows()])
pool.close()
pool.join()

df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SSS_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[((df["7"]>75) & (df["7"]<110)) | ((df["7"]>110) & (df["7"]<165))]
def plot(event,station,time,shift):
    names=["proc_obsd_1D_crust","proc_obsd_3D_crust","proc_obsd_glad","proc_real_data"]
    tags =["proc_obsd_1","proc_obsd_3","proc_obsd_25","proc_real_data"]
    labels=["S40RTS+1D Crust","S40RTS+Crust2.0","GLAD_M25","Real Data"]
    colors=["red","skyblue","green","blueviolet"]
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s.proc_synt.h5")
    data=ds.waveforms[station].proc_synt.select(component="T")[0]
    sampling_rate = data.stats.sampling_rate
    dt            = data.stats.delta
    anltc         = hilbert(data.data)
    envp          = np.abs(anltc)
    arr_ind       = int(time*sampling_rate)
    minus         = int(50*sampling_rate)
    plus          = int(50*sampling_rate)
    p_data        = data.data[arr_ind-minus:arr_ind+plus]
    p_envp        = envp[arr_ind-minus:arr_ind+plus]
    time_series          = data.times()[arr_ind-minus:arr_ind+plus]
    plt.figure(1,figsize=(15,5))
    plt.plot(time_series,p_data,color="black",linewidth=1,label="1D PREM")
    plt.plot(time_series,p_envp,linestyle="--",color="black",linewidth=1,label="1D PREM envelope")
    del ds
    
    for j,n in enumerate(names):
        ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s."+n+".h5")
        data=ds.waveforms[station][tags[j]].select(component="T")[0]
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        anltc         = hilbert(data.data)
        envp          = np.abs(anltc)
        arr_ind       = int((time+shift[j])*sampling_rate)
        minus         = int(50*sampling_rate)
        plus          = int(50*sampling_rate)
        p_data        = data.data[arr_ind-minus:arr_ind+plus]
        p_envp        = envp[arr_ind-minus:arr_ind+plus]
        #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
        del ds
        plt.plot(time_series,p_data,color=colors[j],linewidth=1,label=labels[j]+"(shift = "+str(shift[j])+")")
        plt.plot(time_series,p_envp,linestyle="--",color=colors[j],linewidth=1,label=labels[j]+" envelopes")
    plt.legend(loc=(1.04, 1))
    plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SSS/"+event+"_"+station+".png")
    plt.close()


pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot, [(rows[0],rows[1],rows[8],[rows[9],rows[15],rows[21],rows[27]])  for i,rows in df.iterrows()])
pool.close()
pool.join()

df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/ScS_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[((df["7"]>5) & (df["7"]<25)) | ((df["7"]>50) & (df["7"]<65))]
def plot(event,station,time,shift):
    names=["proc_obsd_1D_crust","proc_obsd_3D_crust","proc_obsd_glad","proc_real_data"]
    tags =["proc_obsd_1","proc_obsd_3","proc_obsd_25","proc_real_data"]
    labels=["S40RTS+1D Crust","S40RTS+Crust2.0","GLAD_M25","Real Data"]
    colors=["red","skyblue","green","blueviolet"]
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s.proc_synt.h5")
    data=ds.waveforms[station].proc_synt.select(component="T")[0]
    sampling_rate = data.stats.sampling_rate
    dt            = data.stats.delta
    anltc         = hilbert(data.data)
    envp          = np.abs(anltc)
    arr_ind       = int(time*sampling_rate)
    minus         = int(50*sampling_rate)
    plus          = int(50*sampling_rate)
    p_data        = data.data[arr_ind-minus:arr_ind+plus]
    p_envp        = envp[arr_ind-minus:arr_ind+plus]
    time_series          = data.times()[arr_ind-minus:arr_ind+plus]
    plt.figure(1,figsize=(15,5))
    plt.plot(time_series,p_data,color="black",linewidth=1,label="1D PREM")
    plt.plot(time_series,p_envp,linestyle="--",color="black",linewidth=1,label="1D PREM envelope")
    del ds
    
    for j,n in enumerate(names):
        ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s."+n+".h5")
        data=ds.waveforms[station][tags[j]].select(component="T")[0]
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        anltc         = hilbert(data.data)
        envp          = np.abs(anltc)
        arr_ind       = int((time+shift[j])*sampling_rate)
        minus         = int(50*sampling_rate)
        plus          = int(50*sampling_rate)
        p_data        = data.data[arr_ind-minus:arr_ind+plus]
        p_envp        = envp[arr_ind-minus:arr_ind+plus]
        #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
        del ds
        plt.plot(time_series,p_data,color=colors[j],linewidth=1,label=labels[j]+"(shift = "+str(shift[j])+")")
        plt.plot(time_series,p_envp,linestyle="--",color=colors[j],linewidth=1,label=labels[j]+" envelopes")
    plt.legend(loc=(1.04, 1))
    plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/ScS/"+event+"_"+station+".png")
    plt.close()


pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot, [(rows[0],rows[1],rows[8],[rows[9],rows[15],rows[21],rows[27]])  for i,rows in df.iterrows()])
pool.close()
pool.join()

df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/Sdiff_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[(df["7"]>100) & (df["7"]<150)]
def plot(event,station,time,shift):
    names=["proc_obsd_1D_crust","proc_obsd_3D_crust","proc_obsd_glad","proc_real_data"]
    tags =["proc_obsd_1","proc_obsd_3","proc_obsd_25","proc_real_data"]
    labels=["S40RTS+1D Crust","S40RTS+Crust2.0","GLAD_M25","Real Data"]
    colors=["red","skyblue","green","blueviolet"]
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s.proc_synt.h5")
    data=ds.waveforms[station].proc_synt.select(component="T")[0]
    sampling_rate = data.stats.sampling_rate
    dt            = data.stats.delta
    anltc         = hilbert(data.data)
    envp          = np.abs(anltc)
    arr_ind       = int(time*sampling_rate)
    minus         = int(50*sampling_rate)
    plus          = int(50*sampling_rate)
    p_data        = data.data[arr_ind-minus:arr_ind+plus]
    p_envp        = envp[arr_ind-minus:arr_ind+plus]
    time_series          = data.times()[arr_ind-minus:arr_ind+plus]
    plt.figure(1,figsize=(15,5))
    plt.plot(time_series,p_data,color="black",linewidth=1,label="1D PREM")
    plt.plot(time_series,p_envp,linestyle="--",color="black",linewidth=1,label="1D PREM envelope")
    del ds
    
    for j,n in enumerate(names):
        ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s."+n+".h5")
        data=ds.waveforms[station][tags[j]].select(component="T")[0]
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        anltc         = hilbert(data.data)
        envp          = np.abs(anltc)
        arr_ind       = int((time+shift[j])*sampling_rate)
        minus         = int(50*sampling_rate)
        plus          = int(50*sampling_rate)
        p_data        = data.data[arr_ind-minus:arr_ind+plus]
        p_envp        = envp[arr_ind-minus:arr_ind+plus]
        #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
        del ds
        plt.plot(time_series,p_data,color=colors[j],linewidth=1,label=labels[j]+"(shift = "+str(shift[j])+")")
        plt.plot(time_series,p_envp,linestyle="--",color=colors[j],linewidth=1,label=labels[j]+" envelopes")
    plt.legend(loc=(1.04, 1))
    plt.savefig("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/Sdiff/"+event+"_"+station+".png")
    plt.close()


pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot, [(rows[0],rows[1],rows[8],[rows[9],rows[15],rows[21],rows[27]])  for i,rows in df.iterrows()])
pool.close()
pool.join()

    
"""