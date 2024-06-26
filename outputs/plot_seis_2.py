import matplotlib.pyplot as plt
import pandas as pd
from pyasdf import ASDFDataSet
from scipy.signal import hilbert
import multiprocessing as mp
import numpy as np
import os
from mpl_toolkits.basemap import Basemap
from scipy import signal

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 22
f_min=1/100.0
f_max=1/40.0

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
file_names      = ["_real","_3D","_glad","_1D"]
ev_list=['200903122323A']

df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/S_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
k=0
df=df[(df["7"]>30) & (df["7"]<70)]

df_plot=df[df[str(16+0*6+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]
for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+2)+df])<12)]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+1)+df])>0.75)]
#df_plot=df_plot[(np.abs(np.log(df_plot[str(16+k*6+5)+file_names[0]]))>1.5)]


def plot(event,station,time,loc):
    lws=[1,1,1,2]
    try:
        names=["proc_obsd_3D_crust","proc_obsd_1D_crust","proc_obsd_glad","proc_real_data"]#,"proc_prem_3D_atten"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
        tags =["proc_obsd_3","proc_obsd_1","proc_obsd_25","proc_real_data"]#"proc_prem_3","proc_obsd_25",
        labels=['S40RTS_3D_crust',"S40RTS_1D_crust","GLAD_M25","real_data"]#"GLAD_M25","PREM_Crust2.0"] 
        colors=["brown","green","cyan","blueviolet"]
        #print("Done")
        ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s.proc_synt.h5",mode="r")
        data=ds.waveforms[station].proc_synt.select(component="T")[0]
        data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        #anltc         = hilbert(data.data)
        #envp          = np.abs(anltc)
        arr_ind       = int(time*sampling_rate)
        minus         = int(500*sampling_rate)
        plus          = int(500*sampling_rate)
        p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        #p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        taper_wdow    = signal.windows.tukey(len(p_data),alpha=0.05)
        p_data        = p_data*taper_wdow
        #p_envp        = p_envp*taper_wdow
        time_series          = data.times()[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        elat          = loc[1]
        elon          = loc[0]
        slat          = loc[4]
        slon         = loc[3]
        EVT           = np.array([0,loc[2]])
        STA           = np.array([loc[5],0])
        
        
        
        fig, (a0, a1) = plt.subplots(2, 1, height_ratios=[2, 1],figsize=(25,15))
        a0.plot(time_series,p_data,color="black",lw=1,label="1D PREM")
        #a0.plot(time_series,p_envp,ls="--",color="black",lw=1)
        del ds
        
        for j,n in enumerate(names):
            ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s."+n+".h5",mode="r")
            data=ds.waveforms[station][tags[j]].select(component="T")[0]
            data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
            sampling_rate = data.stats.sampling_rate
            dt            = data.stats.delta
            #anltc         = hilbert(data.data)
            #envp          = np.abs(anltc)
            arr_ind       = int(time*sampling_rate)
            minus         = int(500*sampling_rate)
            plus          = int(500*sampling_rate)
            p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            #p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            taper_wdow    = signal.windows.tukey(len(p_data),alpha=0.05)
            p_data        = p_data*taper_wdow
            #p_envp        = p_envp*taper_wdow

            time_series          = data.times()[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            del ds
            a0.plot(time_series,p_data,color=colors[j],lw=lws[j],label=labels[j])
            #a0.plot(time_series,p_envp,ls="--",color=colors[j],lw=1)
        a0.legend(loc="upper right")
        a0.arxvline(time)
        map = Basemap(projection='robin', lon_0=-150, resolution="c",ax=a1)
        map.drawmapboundary(fill_color='#cccccc')
        map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)
        map.scatter(elon, elat, s=80, zorder=10,color="red", marker="o",edgecolor="k")
        map.scatter(slon, slat, s=80, zorder=10,color="lime", marker="o",edgecolor="k")
        #if rat>=0:
           #map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='red')
           #plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/S/"+event+"_"+station+".png",dpi=200)
        #else:
        map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='blue')
        plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/S_test/"+event+"_"+station+".png",dpi=200)
        #print(event+" "+station)
        plt.close(fig)
        return 0
    except Exception:
           return 0

print("done")
for i,rows in df_plot.iterrows():
    plot(rows[0],rows[1],rows[8],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]])
#pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
#results=pool.starmap(plot, [(rows[0],rows[1],rows[8],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]])  for i,rows in df_plot.iterrows()])
#pool.close()
#pool.join()

"""
df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SS_all.txt",header=0,delimiter=",",skipinitialspace=True)
k=1
df=df[(df["7"]>50) & (df["7"]<140)]
df_plot=df[df[str(16+k*6+5)+"_real"].notna()]

for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+2)+df])<12)]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+1)+df])>0.75)]
df_plot=df_plot[(np.abs(np.log(df_plot[str(16+k*6+5)+file_names[0]]))>1.5)]

def plot(event,station,time,loc):
    lws=[1,1,1,2]
    try:
        names=["proc_obsd_3D_crust","proc_obsd_1D_crust","proc_obsd_glad","proc_real_data"]#,"proc_prem_3D_atten"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
        tags =["proc_obsd_3","proc_obsd_1","proc_obsd_25","proc_real_data"]#"proc_prem_3","proc_obsd_25",
        labels=['S40RTS_3D_crust',"S40RTS_1D_crust","GLAD_M25","real_data"]#"GLAD_M25","PREM_Crust2.0"] 
        colors=["brown","green","cyan","blueviolet"]
        ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s.proc_synt.h5",mode="r")
        data=ds.waveforms[station].proc_synt.select(component="T")[0]
        data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        #anltc         = hilbert(data.data)
        #envp          = np.abs(anltc)
        arr_ind       = int(time*sampling_rate)
        minus         = int(500*sampling_rate)
        plus          = int(500*sampling_rate)
        p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        #p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        taper_wdow    = signal.windows.tukey(len(p_data),alpha=0.05)
        p_data        = p_data*taper_wdow
       # p_envp        = p_envp*taper_wdow

        time_series          = data.times()[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        elat          = loc[1]
        elon          = loc[0]
        slat          = loc[4]
        slon         = loc[3]
        EVT           = np.array([0,loc[2]])
        STA           = np.array([loc[5],0])

        
        
        fig, (a0, a1) = plt.subplots(2, 1, height_ratios=[2, 1],figsize=(25,15))
        a0.plot(time_series,p_data,color="black",lw=1,label="1D PREM")
        #a0.plot(time_series,p_envp,ls="--",color="black",lw=1)
        del ds
        
        for j,n in enumerate(names):
            ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s."+n+".h5",mode="r")
            data=ds.waveforms[station][tags[j]].select(component="T")[0]
            data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
            sampling_rate = data.stats.sampling_rate
            dt            = data.stats.delta
            #anltc         = hilbert(data.data)
            #envp          = np.abs(anltc)
            arr_ind       = int((time)*sampling_rate)
            minus         = int(500*sampling_rate)
            plus          = int(500*sampling_rate)
            p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            #p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            taper_wdow    = signal.windows.tukey(len(p_data),alpha=0.05)
            p_data        = p_data*taper_wdow
            #p_envp        = p_envp*taper_wdow

            #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
            del ds
            a0.plot(time_series,p_data,color=colors[j],lw=lws[j],label=labels[j])
            #a0.plot(time_series,p_envp,ls="--",color=colors[j],lw=1)
        a0.legend(loc="upper right")
        map = Basemap(projection='robin', lon_0=-150, resolution="c",ax=a1)
        map.drawmapboundary(fill_color='#cccccc')
        map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)
        map.scatter(elon, elat, s=80, zorder=10,color="red", marker="o",edgecolor="k")
        map.scatter(slon, slat, s=80, zorder=10,color="lime", marker="o",edgecolor="k")
        #if rat>=0:
        #   map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='red')
        #   plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/SS2/red/"+event+"_"+station+".png")
        #else:
        map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='blue')
        plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/SS/"+event+"_"+station+".png")
        plt.close(fig)
        return 0
    except ValueError:
           return 0

for i,rows in df_plot.iterrows():
    plot(rows[0],rows[1],rows[8],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]])

df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[((df["7"]>75) & (df["7"]<110)) | ((df["7"]>110) & (df["7"]<165))]
df_plot=df[df[str(16+2*6+5)+"_real"].notna()]
k=2
for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+2)+df])<12)]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+1)+df])>0.75)]
df_plot=df_plot[(np.abs(np.log(df_plot[str(16+k*6+5)+file_names[0]]))>1.5)]
def plot(event,station,time,loc):
    lws=[1,1,1,2]
    try:
        names=["proc_obsd_3D_crust","proc_obsd_1D_crust","proc_obsd_glad","proc_real_data"]#,"proc_prem_3D_atten"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
        tags =["proc_obsd_3","proc_obsd_1","proc_obsd_25","proc_real_data"]#"proc_prem_3","proc_obsd_25",
        labels=['S40RTS_3D_crust',"S40RTS_1D_crust","GLAD_M25","real_data"]#"GLAD_M25","PREM_Crust2.0"] 
        colors=["brown","green","cyan","blueviolet"]
        ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s.proc_synt.h5",mode="r")
        data=ds.waveforms[station].proc_synt.select(component="T")[0]
        data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        #anltc         = hilbert(data.data)
        #envp          = np.abs(anltc)
        arr_ind       = int(time*sampling_rate)
        minus         = int(500*sampling_rate)
        plus          = int(500*sampling_rate)
        p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        #p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        taper_wdow    = signal.windows.tukey(len(p_data),alpha=0.05)
        p_data        = p_data*taper_wdow
        #p_envp        = p_envp*taper_wdow

        time_series          = data.times()[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        elat          = loc[1]
        elon          = loc[0]
        slat          = loc[4]
        slon         = loc[3]
        EVT           = np.array([0,loc[2]])
        STA           = np.array([loc[5],0])

        
        
        fig, (a0, a1) = plt.subplots(2, 1, height_ratios=[2, 1],figsize=(25,15))
        a0.plot(time_series,p_data,color="black",lw=1,label="1D PREM")
        #a0.plot(time_series,p_envp,ls="--",color="black",lw=1)
        del ds
        
        for j,n in enumerate(names):
            ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s."+n+".h5",mode="r")
            data=ds.waveforms[station][tags[j]].select(component="T")[0]
            data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
            sampling_rate = data.stats.sampling_rate
            dt            = data.stats.delta
            #anltc         = hilbert(data.data)
            #envp          = np.abs(anltc)
            arr_ind       = int((time)*sampling_rate)
            minus         = int(500*sampling_rate)
            plus          = int(500*sampling_rate)
            p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            #p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            taper_wdow    = signal.windows.tukey(len(p_data),alpha=0.05)
            p_data        = p_data*taper_wdow
            p_envp        = p_envp*taper_wdow
            
            #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
            del ds
            a0.plot(time_series,p_data,color=colors[j],lw=lws[j],label=labels[j])
            #a0.plot(time_series,p_envp,ls="--",color=colors[j],lw=1)
        a0.legend(loc="upper right")
        map = Basemap(projection='robin', lon_0=-150, resolution="c",ax=a1)
        map.drawmapboundary(fill_color='#cccccc')
        map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)
        map.scatter(elon, elat, s=80, zorder=10,color="red", marker="o",edgecolor="k")
        map.scatter(slon, slat, s=80, zorder=10,color="lime", marker="o",edgecolor="k")
        #if rat>=0:
        #   map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='red')
        #   plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS2/red/"+event+"_"+station+".png")
        map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='blue')
        plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS/"+event+"_"+station+".png")
        plt.close(fig)
        return 0
    except ValueError:
           return 0

for i,rows in df_plot.iterrows():
    plot(rows[0],rows[1],rows[8],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]])


df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/ScS_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[((df["7"]>5) & (df["7"]<25)) | ((df["7"]>50) & (df["7"]<65))]
def plot(event,station,time,shift,loc,rat):
    for j in shift:
        if np.abs(j)>15:
            return None
    try:
        names=["proc_prem_atten_3d","proc_obsd_1D_atten","proc_obsd_3D_atten","proc_real_data"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
        tags =["proc_prem_13","proc_obsd_13","proc_obsd_33","proc_real_data"]#"proc_prem_3","proc_obsd_25",
        labels=["PREM_QRFSI12",'S40RTS_1D_Crust_QRFSI12', 'S40RTS_3D_crust_QRFSI12',"real_data"]
        colors=["brown","cyan","yellowgreen","darkslategrey","blueviolet"]
        ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-040s.proc_synt.h5",mode="r")
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
            ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-040s."+n+".h5",mode="r")
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
           plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/ScS2/red/"+event+"_"+station+".png")
        else:
            map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='blue')
            plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/ScS2/blue/"+event+"_"+station+".png")
        plt.close(fig)
    except ValueError:
           return 0


pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(plot, [(rows[0],rows[1],rows[8],[rows[16],rows[22],rows[28],rows[34],rows[10]],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]],np.log(rows[23])/np.log(rows[17]))  for i,rows in df.iterrows()])
pool.close()
pool.join()

df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/Sdiff_all.txt",header=0,delimiter=",",skipinitialspace=True)
#print(df)
df=df[(df["7"]>100) & (df["7"]<150)]
df_plot=df[df[str(16+4*6+5)+"_real"].notna()]
k=4
for i,df in enumerate(file_names):
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+2)+df])<12)]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+1)+df])>0.75)]
df_plot=df_plot[(np.abs(np.log(df_plot[str(16+k*6+5)+file_names[0]]))>1.5)]
def plot(event,station,time,loc):
    lws=[1,1,1,2]
    try:
        names=["proc_obsd_3D_crust","proc_obsd_1D_crust","proc_obsd_glad","proc_real_data"]#,"proc_prem_3D_atten"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
        tags =["proc_obsd_3","proc_obsd_1","proc_obsd_25","proc_real_data"]#"proc_prem_3","proc_obsd_25",
        labels=['S40RTS_3D_crust',"S40RTS_1D_crust","GLAD_M25","real_data"]#"GLAD_M25","PREM_Crust2.0"] 
        colors=["brown","green","cyan","blueviolet"]
        ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s.proc_synt.h5",mode="r")
        data=ds.waveforms[station].proc_synt.select(component="T")[0]
        data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
        sampling_rate = data.stats.sampling_rate
        dt            = data.stats.delta
        anltc         = hilbert(data.data)
        envp          = np.abs(anltc)
        arr_ind       = int(time*sampling_rate)
        minus         = int(500*sampling_rate)
        plus          = int(500*sampling_rate)
        p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
        p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
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
            ds=ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s."+n+".h5",mode="r")
            data=ds.waveforms[station][tags[j]].select(component="T")[0]
            data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
            sampling_rate = data.stats.sampling_rate
            dt            = data.stats.delta
            #anltc         = hilbert(data.data)
            #envp          = np.abs(anltc)
            arr_ind       = int((time)*sampling_rate)
            minus         = int(500*sampling_rate)
            plus          = int(500*sampling_rate)
            p_data        = data.data[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            #p_envp        = envp[max(0,arr_ind-minus):min(int(3600*sampling_rate),arr_ind+plus)]
            #time_series          = data.times()[arr_ind-minus:arr_ind+plus]
            del ds
            a0.plot(time_series,p_data,color=colors[j],lw=lws[j],label=labels[j])
            #a0.plot(time_series,p_envp,ls="--",color=colors[j],lw=1)
        a0.legend(loc="upper right")
        map = Basemap(projection='robin', lon_0=-150, resolution="c",ax=a1)
        map.drawmapboundary(fill_color='#cccccc')
        map.fillcontinents(color='white', lake_color='#cccccc',zorder=0)
        map.scatter(elon, elat, s=80, zorder=10,color="red", marker="o",edgecolor="k")
        map.scatter(slon, slat, s=80, zorder=10,color="lime", marker="o",edgecolor="k")
        #if rat>=0:
        #   map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='red')
        #   plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/Sdiff2/red/"+event+"_"+station+".png")
        #
        map.drawgreatcircle(elon,elat,slon,slat,ls='None',marker='o',markersize=2,color='blue')
        plt.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/Sdiff/"+event+"_"+station+".png")
        plt.close(fig)
        return 0
    except ValueError:
           return 0

for i,rows in df_plot.iterrows():
    plot(rows[0],rows[1],rows[8],[rows[2],rows[3],rows[4],rows[5],rows[6],rows[7]])

"""