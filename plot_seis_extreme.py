from pyasdf import ASDFDataSet
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
import obspy
import csv
import pygmt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from matplotlib.gridspec import GridSpec


SMALL_SIZE = 24
MEDIUM_SIZE = 24
BIGGER_SIZE = 50
us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"] 

        # controls default text sizes
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


waveforms={}
evlas=[]
evlos=[]
stlas=[]
stlos=[]
rats=[]
events_l=[]
plot_dir="new_figures2"
event_lats=[]
event_lons=[]


file_names=[]
labels=[]
n_list=[872,1482,2405,3326,3941,4534,4922,5405,6026,8112,11289,12181,13332,15821,15916,15988]
for n in n_list:
    labels.append(str(n))
    file_names.append("_"+str(n))
file_names.append("_glad")
file_names.append("_prem")
df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
df = df[~df[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]

# if k==0:
#     df_plot=df_plot[(df_plot["7"]>=40) & (df_plot["7"]<80)]
# elif k==1:
#         df_plot=df_plot[(df_plot["7"]>=80) & (df_plot["7"]<110)]
# elif k==2:
df_plot=df[(df["7"]>=110) & (df["7"]<140)]#df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
for i,df in enumerate(file_names):
    for w in [0,1]:
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]

print(len(df_plot))
df_plot=df_plot[np.abs(np.log(df_plot[str(8)+"_prem"]))>1.0]

df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
print(len(df_plot))
df_plot.to_csv("SS_S_all.txt",index=False)
with open("SS_S_all.txt","r") as f:
    reader=csv.reader(f,delimiter=",",skipinitialspace=True)
    next(reader)
    for row in reader:
        if row[0] not in waveforms.keys():
           waveforms[row[0]]=[row[1]]
        else:
              waveforms[row[0]].append(row[1])
        if row[0] not in events_l:
           events_l.append(row[0])
           event_lats.append(float(row[3]))
           event_lons.append(float(row[2])) 
        evlas.append(float(row[3]))
        evlos.append(float(row[2]))
        stlas.append(float(row[6]))
        stlos.append(float(row[5]))
        rats.append(float(row[8]))

# print(events_l)
# print(event_lats)
# print(event_lons)
fig=pygmt.Figure()

pygmt.makecpt(cmap='polar',reverse=True,series=[-3,3])
fig.basemap(region="d",projection="N-150/12c",frame="ag15")
#fig.coast(shorelines=True, frame=True)   
#fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
for i in range(len(evlas)): 
   fig.plot(x=[evlos[i],stlos[i]],y=[evlas[i],stlas[i]],zvalue=np.log(rats[i]),pen="0.5p,+z",cmap=True)
fig.plot(x=event_lons,y=event_lats,style="c0.2c",fill="green",transparency=40)
fig.coast(shorelines=True, frame=True)
#if i==0:
#   fig.colorbar()
fig.savefig(plot_dir+"/relative/seis/extremes/2/amp_SS_Sdiff.png",dpi=600)
plt.close()
def plot_seis(event,waveforms):
    
    sta_list=waveforms[event]
    tmax=40
    tmin=20
    colors=["gold","red","blue","lime","blue","grey"]

    df_plot=pd.read_csv("SS_S_all.txt",skipinitialspace=True,delimiter=",")
    

    n_list=["1482","8112","15916","GLAD_M25"]
    label =["proc_obsd_3D_crust_1482","proc_obsd_3D_crust","proc_obsd_3D_crust_15916","proc_obsd_glad"]
    tags=["proc_obsd","proc_obsd","proc_obsd","proc_obsd","proc_obsd"]

    # n_list=["8112","PREM","GLAD_M25"]
    # label=["proc_obsd_3D_crust","proc_synt","proc_obsd_glad"]
    # tags=["proc_obsd","proc_synt","proc_obsd"]



    df = ASDFDataSet("/media/ayon/T7/proc/"+event+".T017-250s.proc_real_data.h5",mode="r")
    for origins in df.events[0].origins:
        origin_type=origins.resource_id.id.split("/")[3]
        if origin_type=='origin#cmt':
            second_origin=origins
            break
    evla      = second_origin.latitude
    evlo      = second_origin.longitude
    evdp      = second_origin.depth/1000.0
    for sta in sta_list:
        value = np.log(df_plot[(df_plot["0"]==event) & (df_plot["1"]==sta)][str(8)+"_prem"].values[0])
        #data1=df.waveforms[sta].proc_obsd.select(component="T")[0]
        #data2=df2.waveforms[sta].proc_obsd.select(component="T")[0]
        real=df.waveforms[sta].proc_real_data.select(component="T")[0]
        stax=df.waveforms[sta].StationXML
        #print(stax[0])
        stla=stax[0][0].latitude
        stlo=stax[0][0].longitude

        #print(evla,evlo,stla,stlo)
        print(sta,event)
        st=obspy.Stream()
        st.append(real)
        for d,i in enumerate(label):
            #print(i)
            df2=ASDFDataSet("/media/ayon/T7/proc/"+event+".T017-250s."+i+".h5",mode="r")
            # print(df2.waveforms["PE.IUPA"].proc_obsd.select(component="T"))
            #print(df2.waveforms["CN.FRB"])
            data=df2.waveforms[sta][tags[d]].select(component="T")[0]
            st.append(data)
            del df2
        
       
        st.filter("bandpass",freqmin=1/tmax,freqmax=1/tmin,corners=6,zerophase=True)
        sampling_rate = st[0].stats.sampling_rate

        model=TauPyModel("prem")
        gcarc=locations2degrees(evla,evlo,stla,stlo)
        #print(model.get_travel_times(evdp,gcarc,phase_list=[phase])[0])
        try:
            time=model.get_travel_times(evdp,gcarc,phase_list=["S"])[0].time
        except IndexError:
            time=model.get_travel_times(evdp,gcarc,phase_list=["Sdiff"])[0].time 

        t_bef=1*tmax
        t_aft=(1.25)*tmax
        arr_index=time*sampling_rate
        arr=[arr_index,0]
        start_index  = int(arr_index-t_bef*sampling_rate)
        end_index    = int(arr_index+t_aft*sampling_rate)

        #fig,ax=plt.subplots(1,1,figsize=(9,5))
        #fig1,ax1=plt.subplots(1,1,figsize=(9,5))
        fig=plt.figure(figsize=(15,12))
        gs = GridSpec(3, 2, figure=fig)
        ax=[0,0,0,0]
        ax[2] = fig.add_subplot(gs[0, :])
        ax[0] = fig.add_subplot(gs[1, 0])
        ax[1] = fig.add_subplot(gs[1, 1])
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)
        #ax1.plot(st[0].times(),st[0].data,color="blue",label=n1[1:])
        #ax1.plot(st[1].times(),st[0].data,color="red",label=n2[1:])
        ##ax1.plot(st[2].times(),st[0].data,color="black",label="Real")
        #ax.plot(st[0].times()[start_index:end_index],st[0].data[start_index:end_index],color="blue",label=n1[1:])
        #ax.plot(st[1].times()[start_index:end_index],st[0].data[start_index:end_index],color="red",label=n2[1:])
        data_p,_     = window_taper(st[0][start_index:end_index],1,"cos")
        nrm = max(max(np.abs(stream.data[start_index:end_index])) for stream in st)
        ax[0].plot(st[0].times()[start_index:end_index],data_p/nrm,color="black",label="Observed",lw=2)
        
        for k,n in enumerate(n_list[:-1]):
            data_p,_     = window_taper(st[k+1][start_index:end_index],1,"cos")
            if n!=12181 and n!=4534:
               ax[0].plot(st[0].times()[start_index:end_index],data_p/nrm,color=colors[k],label="N="+str(n),lw=1.5)
            else:
                ax[0].plot(st[0].times()[start_index:end_index],data_p/nrm,color=colors[k],label="N="+str(n),ls="--",lw=1.5)
        data_p,_     = window_taper(st[-1][start_index:end_index],1,"cos")
        ax[0].plot(st[0].times()[start_index:end_index],data_p/nrm,color="grey",label="GLAD-M25",lw=1.5)
        #ax.legend(loc="upper right")
        ax[0].xaxis.set_major_locator(plt.MaxNLocator(5))
        ax[0].grid(True, which='both', linestyle='dotted', linewidth=0.5)
        ax[0].set_title("S")
        ax[0].set_ylim(-1,1)
        ax[0].set_xlim([int(start_index)/4,int(end_index)/4])
        #ax.axis("off")
        ax[0].set_xlabel("Time(s)")
        
        time=model.get_travel_times(evdp,gcarc,phase_list=["SS"])[0].time

        t_bef=1*tmax
        t_aft=(1.25)*tmax
        arr_index=time*sampling_rate
        arr[1]=arr_index
        start_index  = int(arr_index-t_bef*sampling_rate)
        end_index    = int(arr_index+t_aft*sampling_rate)
        # ax2.axvspan(st[0].times()[start_index],st[0].times()[end_index], color='grey', alpha=0.5)
        # ax2.text(st[0].times()[start_index+30],max(st[0].data)/1.3,"SS",fontsize=20)
        nrm = max(max(np.abs(stream.data[start_index:end_index])) for stream in st)
        data_p,_     = window_taper(st[0][start_index:end_index],1,"cos")
        ax[1].plot(st[0].times()[start_index:end_index],data_p/nrm,color="black",label="Observed",lw=2)
        for k,n in enumerate(n_list[:-1]):
            data_p,_     = window_taper(st[k+1][start_index:end_index],1,"cos")
            if n!=12181 and n!=4534:
               ax[1].plot(st[0].times()[start_index:end_index],data_p/nrm,color=colors[k],label="N="+str(n),lw=1.5)
            else:
                ax[1].plot(st[0].times()[start_index:end_index],data_p/nrm,color=colors[k],label="N="+str(n),ls="--",lw=1.5)
        data_p,_     = window_taper(st[-1][start_index:end_index],1,"cos")
        ax[1].plot(st[0].times()[start_index:end_index],data_p/nrm,color="grey",label="GLAD-M25",lw=1.5)
        ax[1].xaxis.set_major_locator(plt.MaxNLocator(5))
        ax[1].set_ylim(-1,1)
        ax[1].grid(True, which='both', linestyle='dotted', linewidth=1)
        #ax1.legend(loc="upper right")
        ax[1].set_title("SS")
        ax[1].set_xlabel("Time(s)")
        ax[1].set_xlim([int(start_index)/4,int(end_index)/4])
        #ax1.axis("off")
        plt.rc('xtick', labelsize=24)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=24)
        start_index  = int(arr[0]-200*sampling_rate)
        end_index    = int(arr[1]+800*sampling_rate)
        
        lims         =  max(max(np.abs(stream.data[start_index:end_index])) for stream in st)
        print(lims)
        st=obspy.Stream()
        st.append(real)
        for d,i in enumerate(label):
            #print(i)
            df2=ASDFDataSet("/media/ayon/T7/proc/"+event+".T017-250s."+i+".h5",mode="r")
            # print(df2.waveforms["PE.IUPA"].proc_obsd.select(component="T"))
            #print(df2.waveforms["CN.FRB"])
            data=df2.waveforms[sta][tags[d]].select(component="T")[0]
            st.append(data)
            del df2
        st.filter("bandpass",freqmin=1/tmax,freqmax=1/tmin,corners=6,zerophase=True)
        sampling_rate = st[0].stats.sampling_rate
        data_p     = st[0][start_index:end_index]
        ax[2].plot(st[0].times()[start_index:end_index],data_p,color="black",label="Observed",lw=1.5)
        for k,n in enumerate(n_list[:-1]):
            data_p     = st[k+1][start_index:end_index]
            if n!=12181 and n!=4534:
               ax[2].plot(st[0].times()[start_index:end_index],data_p,color=colors[k],label="N="+str(n),lw=1)
            else:
                ax[2].plot(st[0].times()[start_index:end_index],data_p,color=colors[k],label="N="+str(n),ls="--",lw=1)
        data_p     = st[-1][start_index:end_index]
        ax[2].plot(st[0].times()[start_index:end_index],data_p,color="grey",label="GLAD-M25",lw=1)
        ax[2].legend(loc="lower right", fontsize=14)
        ax[2].axvspan(st[0].times()[int(arr[0]-t_bef*sampling_rate)],st[0].times()[int(arr[0]+t_aft*sampling_rate)], color='grey', alpha=0.5)
        ax[2].text(st[0].times()[int(arr[0]-t_bef*sampling_rate)],max(st[0].data[start_index:end_index])*0.8,"Sdiff",fontsize=28)
        ax[2].axvspan(st[0].times()[int(arr[1]-t_bef*sampling_rate)],st[0].times()[int(arr[1]+t_aft*sampling_rate)], color='grey', alpha=0.5)
        ax[2].text(st[0].times()[int(arr[1]+10)],max(st[0].data[start_index:end_index])*0.8,"SS",fontsize=28)
        ax[2].set_title("ln(R)="+str(np.round(value,2))+" Depth="+str(evdp)+" km"+" Dist="+str(np.round(gcarc,0)),fontsize=32)
        ax[2].set_xlabel("Time(s)", fontsize=32)
        ax[2].set_ylabel("Displacement(nm)",fontsize=32)
        ax[2].set_ylim(-1.1*lims,1.1*lims)
        #ax2.xaxis.set_major_locator(plt.MaxNLocator(7))
        ax[2].grid(True, which='both', linestyle='dotted', linewidth=1)
        ax[2].set_xlim([int(start_index)/4,int(end_index)/4])
        #ax[2].set_xlim([1000,2500])
        ax_c= fig.add_subplot(gs[2, :], projection=ccrs.PlateCarree(central_longitude=-150))
        ax_c.set_global()
        ax_c.coastlines()
        ax_c.add_feature(cfeature.LAND)
        # ax_c.scatter(evlo, evla, 200, marker="*", color="red", transform=ccrs.PlateCarree(), zorder=100)
        ax_c.scatter(evlo, evla, 200, marker="*", color="red", transform=ccrs.PlateCarree(), zorder=100)
        ax_c.plot([evlo, stlo], [evla, stla], "k", transform=ccrs.Geodetic())
        

        #fig.tight_layout()
        #fig1.tight_layout()
        fig.tight_layout()
        plot_dir="new_figures2"
        fig.savefig(plot_dir+"/relative/seis/extremes/2/"+str(np.round(value,2))+"SS_S"+"_"+event+"_"+sta+"_"+str(np.round(gcarc,1))+".pdf")
        #fig1.savefig(plot_dir+"/relative/seis/north_pacific/2/SS"+"_"+event+"_"+sta+"_"+str(np.round(gcarc,1))+".pdf", dpi=600)
        #fig.savefig(plot_dir+"/relative/seis/north_pacific/2/S"+"_"+event+"_"+sta+"_"+str(np.round(gcarc,1))+".pdf",dpi=600)
        #plt.close(fig)
        plt.close(fig)

        st=obspy.Stream()
        real=df.waveforms[sta].proc_real_data.select(component="T")[0]
        st.append(real)
        for d,i in enumerate(label):
            #print(i)
            df2=ASDFDataSet("/media/ayon/T7/proc/"+event+".T017-250s."+i+".h5",mode="r")
            # print(df2.waveforms["PE.IUPA"].proc_obsd.select(component="T"))
            #print(df2.waveforms["CN.FRB"])
            data=df2.waveforms[sta][tags[d]].select(component="T")[0]
            st.append(data)
            del df2
        st_Z=obspy.Stream()
        real=df.waveforms[sta].proc_real_data.select(component="Z")[0]
        st_Z.append(real)
        for d,i in enumerate(label):
            #print(i)
            df2=ASDFDataSet("/media/ayon/T7/proc/"+event+".T017-250s."+i+".h5",mode="r")
            # print(df2.waveforms["PE.IUPA"].proc_obsd.select(component="T"))
            #print(df2.waveforms["CN.FRB"])
            data=df2.waveforms[sta][tags[d]].select(component="Z")[0]
            st_Z.append(data)
            del df2
        st_R=obspy.Stream()
        real=df.waveforms[sta].proc_real_data.select(component="R")[0]
        st_R.append(real)
        for d,i in enumerate(label):
            #print(i)
            df2=ASDFDataSet("/media/ayon/T7/proc/"+event+".T017-250s."+i+".h5",mode="r")
            # print(df2.waveforms["PE.IUPA"].proc_obsd.select(component="T"))
            #print(df2.waveforms["CN.FRB"])
            data=df2.waveforms[sta][tags[d]].select(component="R")[0]
            st_R.append(data)
            del df2
        time=model.get_travel_times(evdp,gcarc,phase_list=["Sdiff"])[0].time
        time2=model.get_travel_times(evdp,gcarc,phase_list=["SS"])[0].time
        fig,ax=plt.subplots(3,1,figsize=(15,13))
        for i,sts in enumerate([st,st_Z,st_R]):
            sts.filter("bandpass",freqmin=1/tmax,freqmax=1/tmin,corners=6,zerophase=False)
            # for tr in sts:
            ax[i].plot(sts[0].times()[2000:10000],sts[0].data[2000:10000],color="black",label="Observed")
            ax[i].plot(sts[-1].times()[2000:10000],sts[1].data[2000:10000],color="red",label="GLAD-M25")
        ax[0].set_title("T")
        ax[1].set_title("Z")
        ax[2].set_title("R")
        ax[0].axvline(x=time, color='blue',lw=1)
        ax[1].axvline(x=time, color='blue',lw=1)
        ax[2].axvline(x=time, color='blue',lw=1)
        ax[0].axvline(x=time2, color='blue',lw=1)
        ax[1].axvline(x=time2, color='blue',lw=1)
        ax[2].axvline(x=time2, color='blue',lw=1)
        plt.savefig("new_figures2/relative/seis/extremes/2/"+str(np.round(value,2))+"SS_S"+"_"+event+"_"+sta+"_"+str(np.round(gcarc,1))+"_ZRT.pdf")
        plt.close(fig)
    del df

for keys in waveforms.keys():
    plot_seis(keys,waveforms)