import matplotlib.pyplot as plt
import pandas as pd
from pyasdf import ASDFDataSet
from scipy.signal import hilbert
from multiprocessing import Pool
import numpy as np
import os
import csv
from mpl_toolkits.basemap import Basemap
from scipy import signal
from matplotlib.lines import Line2D
from obspy.taup import TauPyModel

plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/today"
def window_taper(signal, taper_percentage, taper_type):
    npts = len(signal)
    frac = int(npts * taper_percentage / 2.0 + 0.5)

    idx1 = frac
    idx2 = npts - frac

    if taper_type == "hann":
        window = 0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(0, 2 * frac) / (2 * frac - 1))
    elif taper_type == "cos":
        window = np.cos(np.pi * np.arange(0, 2 * frac) / (2 * frac - 1) - np.pi / 2.0)
    elif taper_type == "cos_p10":
        window = 1. - np.cos(np.pi * np.arange(0, 2 * frac) / (2 * frac - 1)) ** 10

    signal[:idx1] *= window[:frac]
    signal[idx2:] *= window[frac:]

    # return test signal
    test = np.ones(npts)
    test[:idx1] = window[:frac]
    test[idx2:] = window[frac:]

    return signal, test

SMALL_SIZE = 24
MEDIUM_SIZE = 32
BIGGER_SIZE = 38
f_min=1/100.0
f_max=1/40.0

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


ev_list=["201705091352A","201211122042A","201801310707A","201710101853A"]

with open("./../event_list","r") as txt:
     reader=csv.reader(txt, skipinitialspace=True,delimiter=" ")
     for row in reader:
         ev_list.append((row[0]))
     
file_names=["obsd_3D_crust","obsd_1D_crust","obsd_glad","prem_3D_crust","synt"]


# df=[]
# for i,f in enumerate(file_names):
#     temp    = pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/0_"+f+"_real.txt",header=None,delimiter=" ",skipinitialspace=True)
#     df.append(temp)


# d_all=pd.merge(df[1],df[0],how='left',on=[0,1,2,3,4,5,6,7,18,29,40,62],suffixes=("_1D","_3D"))
# d_all=pd.merge(d_all,df[2],how='left',on=[0,1,2,3,4,5,6,7,18,29,40,62])
# d_all=pd.merge(d_all,df[3],how='left',on=[0,1,2,3,4,5,6,7,18,29,40,62],suffixes=("_glad","_prem_3D"))
# d_all=pd.merge(d_all,df[4],how='left',on=[0,1,2,3,4,5,6,7,18,29,40,62])
# d_all=pd.merge(d_all,df[4],how='left',on=[0,1,2,3,4,5,6,7,18,29,40,62],suffixes=("_synt","_fow"))

# names2=["_synt","_3D","_glad"]
# col_list=[]
# for i in range(8):
#     col_list.append(i)
# n_list=[62]
# for j in n_list:
#     col_list.append(j)

# k_l=[4]
# for k in k_l:
    
#     for n in names2:
#         col_list.append(str(18+k*11+2)+n)
#         col_list.append(str(18+k*11+6)+n)
# k=[4]
df_plot =pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/Sdiff_all_real.txt",delimiter=",",skipinitialspace=True)
lws=[1,1,1,1,1,2]
names=["proc_synt","proc_obsd_3D_crust","proc_obsd_glad","proc_real_data"]#,"proc_prem_3D_atten"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
tags =["proc_synt","proc_obsd_1","proc_obsd_3","proc_obsd_25","proc_prem_3","proc_real_data"]#"proc_prem_3","proc_obsd_25",
labels=["real_data","S40RTS_Crust2.0","GLAD_M25"]#"GLAD_M25","PREM_Crust2.0"] 
colors=['black','red','deepskyblue','green','forestgreen','indianred','grey','peru']
colors=["black","red","red","royalblue","forestgreen","red"]
phase=['Sdiff']

names=["proc_real_data","proc_obsd_3D_crust","proc_obsd_glad"]#,"proc_prem_3D_atten"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
tags =["proc_real_data","proc_obsd_3","proc_obsd_25"]#
names2=["_real","_3D","_glad"]

def plot(event,df_plot,station,t_min,t_max,ph):
    df_event=df_plot[(df_plot['0']==event) & (df_plot['1']==station)]
    print(df_event)
    model=TauPyModel(model="prem")
    f_min = 1/t_max
    f_max = 1/t_min
    t_bef = 0.5/f_min
    t_aft = (1+0.5/2)/f_min+100
    
        #print("Done")
        
        #dss=[ds_synt,ds_1D,ds_3D,ds_glad,ds_prem_3D,ds_real]
            
    start_time=1000000
    end_time=0
    fig=plt.figure(figsize=(30,10))
    for l,d in enumerate(names[0:2]):
        
        ds_synt       = ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s."+d+".h5",mode="r")
        
        
        if station in ds_synt.waveforms.list() and len(ds_synt.waveforms[station][tags[l]].select(component="T"))>0:
            # evla=ds_synt.events[0].origins[0].latitude
            # evlo=ds_synt.events[0].origins[0].longitude
            # evdp=ds_synt.events[0].origins[0].depth/1000
            # stla=ds_synt.waveforms[station].StationXML[0][0].latitude
            # stlo=ds_synt.waveforms[station].StationXML[0][0].longitude
            try:
                evla=df_event['3'].iloc[0]
                evlo=df_event['2'].iloc[0]
                evdp=df_event['4'].iloc[0]
                stla=df_event['6'].iloc[0]
                stlo=df_event['5'].iloc[0]
            except IndexError:
                evla=ds_synt.events[0].origins[0].latitude
                evlo=ds_synt.events[0].origins[0].longitude
                evdp=ds_synt.events[0].origins[0].depth/1000
                stla=ds_synt.waveforms[station].StationXML[0][0].latitude
                stlo=ds_synt.waveforms[station].StationXML[0][0].longitude
        
            
            arrival_time=model.get_travel_times_geo(evdp,evla,evlo,stla,stlo,phase_list=['Sdiff'])[0].time
            arrival_time2=model.get_travel_times_geo(evdp,evla,evlo,stla,stlo,phase_list=['SSS'])[0].time
            arrival_time3=model.get_travel_times_geo(evdp,evla,evlo,stla,stlo,phase_list=['SS'])[0].time
            #except IndexError:
            #       print("Phase Not Available")
            #       return 0
            data          =           ds_synt.waveforms[station][tags[l]].select(component="T")[0]
            data          = data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=6,zerophase=True)
            sampling_rate = data.stats.sampling_rate
            dt            = data.stats.delta
            if l==0:
                arr_index       = int(arrival_time*sampling_rate)
                arr_index2       = int(arrival_time2*sampling_rate)
                start_index  = int(arr_index-t_bef*sampling_rate)
                end_index    = int(arr_index2+t_aft*sampling_rate)
                data_wdow    = data.data[start_index:end_index]
                data_wdow,_  = window_taper(data_wdow, 1, "cos")
                time_series  = data.times()[start_index:end_index]
                ax1=fig.add_subplot(1,1,1) 
                        #print("ax1 created")
                ax1.plot(time_series,data_wdow,color=colors[l], lw=2, label=labels[l])
                #rms_synt     = np.sqrt(np.mean(data_wdow**2))
                
            if l!=0:
                #if  row[str(18+k[i]*11+2)+names2[l-1]]==row[str(18+k[i]*11+2)+names2[l-1]]:
                if 2==2:
                    #print(l)  
                    
                    #arr_index     = int((row[j]-row[str(18+k[i]*11+2)+names2[l-1]])*sampling_rate)
                    arr_index       = int(arrival_time*sampling_rate)
                    arr_index2       = int(arrival_time2*sampling_rate)
                    
                    start_index  = int(arr_index-t_bef*sampling_rate)
                    end_index    = int(arr_index2+t_aft*sampling_rate)
                    data_wdow    = data.data[start_index:end_index]
                    data_wdow,_  = window_taper(data_wdow, 1, "cos")
                    time_series  = data.times()[start_index:end_index]
                    #rms_obsd     = np.sqrt(np.mean(data_wdow**2))
                    #print(np.log(rms_obsd/rms_synt))
                        
                    
                    #ax1.plot(time_series,data_wdow,color=colors[l], lw=lws[l], label=str(np.round(df_event['68'+str(names2[l])].iloc[0],2)))                        
                    ax1.plot(time_series,data_wdow,color=colors[l], lw=3, label=labels[l]) 
            ax1.grid(True)
            #ax1.text(arrival_time, 0.5*(max(data_wdow)), 'Sdiff', fontsize=12, color='red')      
            #ax1.text(arrival_time2, 0.5*(max(data_wdow)), 'SS', fontsize=12, color='red') 
            #ax1.text(arrival_time3, 0.5*(max(data_wdow)), 'SSS', fontsize=12, color='red')    
            #ax1.title.set_text(ph)
            ax1.set_xlabel('Time(s)')
            ax1.set_ylabel('Displacement(mm)')
            
            #ax1.set_ylim(-9e-7,9e-7)
            ax1.legend(loc="upper right",frameon=False)
    
            del ds_synt
                    
                
    
    
    # handles=[]
    # for x,df in enumerate(names):
    #     if x==3:
    #         handles.append(Line2D([0], [0], color=colors[x], lw=lws[x], label=labels[x]))
    #     else:
        
    #         handles.append(Line2D([0], [0], color=colors[x], lw=lws[x], label=labels[x]+" "+str(df_event['68'+str(names2[x])])))
    # ax2=fig.add_subplot(2,1,2)
    # ax2.set_axis_off()
    # ax2.legend(handles=handles, loc='center')
    ax1.text(arrival_time, 0.2*(max(data_wdow)), 'Sdiff', fontsize=36, color='blue')
    ax1.text(arrival_time3, 0.2*(max(data_wdow)), 'SS', fontsize=36, color='blue')
    ax1.text(arrival_time2, 0.2*(max(data_wdow)), 'SSS', fontsize=36, color='blue')
    #ax1.set_xlim(
    plt.tight_layout()
    #plt.show()
    plt.savefig(plot_dir+"/"+event+"_"+station+"_"+ph+".png")
    plt.close()
    #print("Saved")
    
    return 0
events=["201906281551A","201705091352A","201704031740A"]
stations=[["LD.SDMD","IU.MACI","LD.FOR"],["GT.LBTB","PR.HUMP","WI.MPOM"],["AK.PAX","N4.R50A","BK.HUMO","NZ.ODZ"]]
events=["201705091352A"]
stations=[["NO.SPA0"]]

for w,event in enumerate(events):
    for station in stations[w]:    
        
        
        t_min=20.0
        t_max=50.0
        ph='Sdiff'
        plot(event,df_plot,station,t_min,t_max,'SS')
    #plot(event,df_plot,station,t_min,t_max,'SS')
    #plot(event,df_plot,station,t_min,t_max,'SSS')

handles=[]
for i,df in enumerate(names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lws[i], label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends")
plt.close()




            
            

