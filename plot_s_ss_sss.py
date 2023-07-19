#importing necessary modules
import os
from obspy.taup import TauPyModel
from obspy import read
import matplotlib.pyplot as plt
import csv
from obspy import Stream 
import numpy as np
import operator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import signal
from scipy import fft
import pygmt
import scipy
import glob
import multiprocessing as mp
import matplotlib.pyplot as plt


def spec_amp(signal,f,delta):
    
    sig_fft = fft.rfft(signal)

    amp_sig=np.real(sig_fft)
    freq=fft.rfftfreq(len(signal),d=delta)
    freq_one=freq[:int(len(sig_fft)/2)]
    disc_freq=1/(len(signal)*delta)
    i=int(f/disc_freq)
    a=0
    k=0
    while freq[i]<=0.05:
          a+=np.abs(amp_sig[i])
          k+=1
          i+=1
    #print(k)
    a=a/float(k)
    return a
    
    
fr=0.03
#Reading the stations where we have data

model=TauPyModel(model='prem')
os.chdir('/home/ayon/WORK/new_events/synthetics/1D_shallow/')

def processing(run):
    print("Working on "+run+" directory")
    dir_1d='/home/ayon/WORK/new_events/synthetics/1D_shallow/'+run+'/OUTPUT_FILES'
    dir_3d='/home/ayon/WORK/new_events/synthetics/3D_crust_shaloow/'+run+'/OUTPUT_FILES'
    stations=[]
    distance=[]
    with open(dir_1d+'/output_list_stations.txt','r') as txt:
         for i in range(3):
             next(txt)
         reader=csv.reader(txt,skipinitialspace=True,delimiter=' ')
         for row in reader:
             stations.append(row[0])
             distance.append(float(row[3]))
    try:
        os.mkdir(dir_3d+'/../figures')        
    except FileExistsError: 
           print("File already exists")
    


    try:  
        with open(dir_3d+'/../stations_done.txt','r') as txt:
             last_line=txt.readlines()[-1]
             start=last_line.split(" ")[0]
    except FileNotFoundError:
           start=0
   #print(start,run)
    for sta in range(int(start),len(stations)):
        
        st=Stream()
        
        try:
            st=read(dir_1d+'/'+str(stations[sta])+'.MXT.sem.sac')
        except FileNotFoundError:
               continue
        try:
            st+=read(dir_3d+'/'+str(stations[sta])+'.MXT.sem.sac')
        except FileNotFoundError:
               continue
        #print(st)
        if distance[sta]<180.0 and distance[sta]>30.0:

           st.detrend()
           st.filter("bandpass",freqmin=1/55.0,freqmax=1/17.0,zerophase=True)
          
        #st.plot()

           evdp=st[0].stats.sac.evdp
           stdp=st[0].stats.sac.stdp
           evla=st[0].stats.sac.evla
           evlo=st[0].stats.sac.evlo
           stla=st[0].stats.sac.stla
           stlo=st[0].stats.sac.stlo
           gcarc=st[0].stats.sac.gcarc
           dt=st[0].stats.sac.delta
           dt3=st[1].stats.sac.delta
           starttime=st[0].stats.starttime
           baz=st[0].stats.sac.baz
        
           freq=4
           delta=1.0/freq
           npts=5400*freq
           for tr in st:
               tr.interpolate(sampling_rate=freq,starttime=starttime,npts=npts)

        
           arrivals=model.get_travel_times(source_depth_in_km=evdp,receiver_depth_in_km=stdp,distance_in_degree=gcarc,phase_list=['S','SS','SSS','Sdiff','ScS'])
           a=len(arrivals)
           #print(arrivals)
##Putting arrivals in the traces header
           s_arrival=0
           s_count=0
           ss_arrival=0
           ss_count=0
           sss_arrival=0
           sss_count=0
           scs_arrival=0
           scs_count=0
           sdiff_arrival=0
           sdiff_count=0
           for j in range(a):
               if arrivals[j].name=='S':
                  if s_arrival==0:
                     s_arrival=arrivals[j].time
                   #print(s_arrival)
                     s_count+=1
               elif arrivals[j].name=='SS':
                    if ss_arrival==0:
                       ss_arrival=arrivals[j].time
                     #print(ss_arrival)
                       ss_count+=1
               elif arrivals[j].name=='SSS':
                    if sss_arrival==0:
                       sss_arrival=arrivals[j].time
                     #print(sss_arrival)
                       sss_count+=1
               elif arrivals[j].name=='ScS':
                    if scs_arrival==0:
                       scs_arrival=arrivals[j].time
                       scs_count+=1
               elif arrivals[j].name=='Sdiff':
                   if sdiff_arrival==0:
                       sdiff_arrival=arrivals[j].time
                       sdiff_count+=1
           for i in range(2):
               if s_count==1:
                  st[i].stats.sac.t0=s_arrival
               # print(s_arrival)
               if ss_count==1:
                  st[i].stats.sac.t1=ss_arrival
               if sss_count==1:
                  st[i].stats.sac.t2=sss_arrival
               if scs_count==1:
                  st[i].stats.sac.t3=scs_arrival
               if sdiff_count==1:
                  st[i].stats.sac.t4=sdiff_arrival



     

##Cross-correlating S, SS and SSS traces
           data_1d=st[0].data
           data_3d=st[1].data
           time_1d=st[0].times()
           time_3d=st[1].times()
           #print(arrivals)
           s_window_start=0
           s_window_end=0
           ss_window_start=0
           scs_window_start=0
           sdiff_window_start=0
           ss_window_end=0
           sss_window_start=0
           sss_window_end=0
           scs_window_end=0
           sss_window_end=0
           synt_anltc = signal.hilbert(data_1d)
           obsd_anltc = signal.hilbert(data_3d)
           obsd_envp = np.abs(obsd_anltc)
           synt_envp = np.abs(synt_anltc)

           obsd_phse = np.unwrap(np.angle(obsd_anltc))
           synt_phse = np.unwrap(np.angle(synt_anltc))


## Doing Calculations
           try:
           ## cutting windows and tapering and calculating shift
               s_window_start=st[1].stats.sac.t0-50.0
               s_window_end=s_window_start+100.0
               index_s_start=int(s_window_start/delta)
               index_s_end=int(s_window_end/delta)
               s_w_1d=data_1d[index_s_start:index_s_end]
               taper_1d=signal.tukey(len(s_w_1d),alpha=0.5)
               s_envp_1d = synt_envp[index_s_start:index_s_end]*taper_1d
               s_w_1d=s_w_1d*taper_1d
               s_w_3d=data_3d[index_s_start:index_s_end]
               taper_3d=signal.tukey(len(s_w_3d),alpha=0.5)
               s_w_3d=s_w_3d*taper_3d
               s_corr=signal.correlate(s_w_1d,s_w_3d,method='auto',mode='full')
               max_ind=np.argmax(s_corr)
               s_shift=max_ind-len(s_w_1d)+1
               shifted_3d_s=data_3d[index_s_start-s_shift:index_s_end-s_shift]
               s_envp_3d = obsd_envp[index_s_start-s_shift:index_s_end-s_shift]*taper_3d
               shifted_3d_s=taper_3d*shifted_3d_s
               index_1d_s,s_max_1d=max(enumerate(np.abs(s_w_1d)),key=operator.itemgetter(1))
               index_3d_s,s_max_3d=max(enumerate(np.abs(shifted_3d_s)),key=operator.itemgetter(1))
            #Calculating Energy and other Statistics
               energy_s=scipy.integrate.trapezoid((s_w_1d-shifted_3d_s)**2,dx=delta)
               s1d_sqr=scipy.integrate.trapezoid(s_w_1d*s_w_1d,dx=delta)
               s3d_sqr=scipy.integrate.trapezoid(shifted_3d_s*shifted_3d_s,dx=delta)
               s1d_s3d=scipy.integrate.trapezoid(s_w_1d*shifted_3d_s,dx=delta)
               s_3d=scipy.integrate.trapezoid(shifted_3d_s**2,dx=delta)
               norm_s=scipy.integrate.trapezoid(s_w_1d**2,dx=delta)
               A_s=energy_s/norm_s
               R_s=(1/len(s_w_1d))*(s_3d/norm_s)
           #Calculating Spectral Amplitudes in 3D and 1D
               
              # for fr in freqlist:
               amp_1d_s=spec_amp(s_w_1d,fr,delta)
               
               amp_3d_s=spec_amp(shifted_3d_s,fr,delta)

               envp_mf_s = np.log(max(s_envp_3d)/max(s_envp_1d))
               max_s_envp_3d=np.argmax(s_envp_3d)
               max_s_envp_1d=np.argmax(s_envp_1d)
               #envp_mf_s = 0.5*scipy.integrate.trapezoid(lne**2)
                         
           except AttributeError:
                  print("No S Phase detected")
           try:
               ss_window_start=st[1].stats.sac.t1-50.0
               ss_window_end=ss_window_start+100.0
               index_ss_start=int(ss_window_start/delta)
               index_ss_end=int(ss_window_end/delta)
               ss_w_1d=data_1d[index_ss_start:index_ss_end]
               ss_envp_1d=synt_envp[index_ss_start:index_ss_end]*signal.tukey(len(ss_w_1d),alpha=0.5)
               ss_w_1d=ss_w_1d*signal.tukey(len(ss_w_1d),alpha=0.5)
               ss_w_3d=data_3d[index_ss_start:index_ss_end]
               ss_w_3d=ss_w_3d*signal.tukey(len(ss_w_3d),alpha=0.5)
               ss_corr=signal.correlate(ss_w_1d,ss_w_3d,method='auto',mode='full')
               max_ind=np.argmax(ss_corr)
               ss_shift=max_ind-len(ss_w_1d)+1
               shifted_3d_ss=data_3d[index_ss_start-ss_shift:index_ss_end-ss_shift]*signal.tukey(len(ss_w_3d),alpha=0.5)
               ss_envp_3d = obsd_envp[index_ss_start-ss_shift:index_ss_end-ss_shift]*signal.tukey(len(ss_w_3d),alpha=0.5)
               index_1d_ss,ss_max_1d=max(enumerate(np.abs(ss_w_1d)),key=operator.itemgetter(1))
               index_3d_ss,ss_max_3d=max(enumerate(np.abs(shifted_3d_ss)),key=operator.itemgetter(1))
       #     for i in range(len(ss_w_1d)):
               energy_ss=scipy.integrate.trapezoid((ss_w_1d-shifted_3d_ss)**2,dx=delta)
               ss1d_sqr=scipy.integrate.trapezoid(ss_w_1d*ss_w_1d,dx=delta)
               ss3d_sqr=scipy.integrate.trapezoid(shifted_3d_ss*shifted_3d_ss,dx=delta)
               ss1d_ss3d=scipy.integrate.trapezoid(ss_w_1d*shifted_3d_ss,dx=delta)
               ss_3d=scipy.integrate.trapezoid(shifted_3d_ss**2,dx=delta)
               norm_ss=scipy.integrate.trapezoid(ss_w_1d**2,dx=delta)
               A_ss=energy_ss/norm_ss
               R_ss=(1/len(ss_w_1d))*(ss_3d/norm_ss)
           #Calculating Spectral Amplitudes in 3D and 1D
               amp_1d_ss=spec_amp(ss_w_1d,fr,delta)
               
               amp_3d_ss=spec_amp(shifted_3d_ss,fr,delta)

               envp_mf_ss = np.log(max(ss_envp_3d)/max(ss_envp_1d))
               max_ss_envp_3d=np.argmax(ss_envp_3d)
               max_ss_envp_1d=np.argmax(ss_envp_1d)

               #envp_mf_ss = 0.5*scipy.integrate.trapezoid(lne**2)

           except AttributeError:
                  print("No SS Phase detected")
           try:
               sss_window_start=st[1].stats.sac.t2-50.0
               sss_window_end=sss_window_start+100.0
               index_sss_start=int(sss_window_start/delta)
               index_sss_end=int(sss_window_end/delta)
               sss_w_1d=data_1d[index_sss_start:index_sss_end]
               sss_envp_1d=synt_envp[index_sss_start:index_sss_end]
               sss_w_1d=sss_w_1d*signal.tukey(len(sss_w_1d),alpha=0.5)            
               sss_w_3d=data_3d[index_sss_start:index_sss_end]*signal.tukey(len(sss_w_1d),alpha=0.5)
               sss_w_3d=sss_w_3d*signal.tukey(len(sss_w_3d),alpha=0.5)
               sss_corr=signal.correlate(sss_w_1d,sss_w_3d,method='auto',mode='full')
               max_ind=np.argmax(sss_corr)
               sss_shift=max_ind-len(sss_w_1d)+1
               shifted_3d_sss=data_3d[index_sss_start-sss_shift:index_sss_end-sss_shift]*signal.tukey(len(sss_w_3d),alpha=0.5)
               sss_envp_3d = obsd_envp[index_sss_start-sss_shift:index_sss_end-sss_shift]*signal.tukey(len(sss_w_3d),alpha=0.5)
               index_1d_sss,sss_max_1d=max(enumerate(np.abs(sss_w_1d)),key=operator.itemgetter(1))
               index_3d_sss,sss_max_3d=max(enumerate(np.abs(shifted_3d_sss)),key=operator.itemgetter(1))
           # for i in range(len(sss_w_1d)):
               energy_sss=scipy.integrate.trapezoid((sss_w_1d-shifted_3d_sss)**2,dx=delta)
               sss1d_sqr=scipy.integrate.trapezoid(sss_w_1d*sss_w_1d,dx=delta)
               sss3d_sqr=scipy.integrate.trapezoid(shifted_3d_sss*shifted_3d_sss,dx=delta)
               sss1d_sss3d=scipy.integrate.trapezoid(sss_w_1d*shifted_3d_sss,dx=delta)
               sss_3d=scipy.integrate.trapezoid(shifted_3d_sss**2,dx=delta)
               norm_sss=scipy.integrate.trapezoid(sss_w_1d**2,dx=delta)
               A_sss=energy_sss/norm_sss
               R_sss=(1/len(sss_w_1d))*(sss_3d/norm_sss)
               amp_1d_sss=spec_amp(sss_w_1d,fr,delta)
               
               amp_3d_sss=spec_amp(shifted_3d_sss,fr,delta)               
               ### Envelope Misfit for SSS
               

               envp_mf_sss = np.log(max(sss_envp_3d)/max(sss_envp_1d))
               max_sss_envp_3d=np.argmax(sss_envp_3d)
               max_sss_envp_1d=np.argmax(sss_envp_1d)

               #envp_mf_sss = 0.5*scipy.integrate.trapezoid(lne**2)
               
               
           except AttributeError:
                  print("No SSS Phase detected")

           try:
               scs_window_start=st[1].stats.sac.t3-50.0
               scs_window_end=scs_window_start+100.0
               index_scs_start=int(scs_window_start/delta)
               index_scs_end=int(scs_window_end/delta)
               scs_w_1d=data_1d[index_scs_start:index_scs_end]
               scs_envp_1d=synt_envp[index_scs_start:index_scs_end]
               scs_w_1d=scs_w_1d*signal.tukey(len(scs_w_1d),alpha=0.5)            
               scs_w_3d=data_3d[index_scs_start:index_scs_end]*signal.tukey(len(scs_w_1d),alpha=0.5)
               scs_w_3d=scs_w_3d*signal.tukey(len(scs_w_3d),alpha=0.5)
               scs_corr=signal.correlate(scs_w_1d,scs_w_3d,method='auto',mode='full')
               max_ind=np.argmax(scs_corr)
               scs_shift=max_ind-len(scs_w_1d)+1
               shifted_3d_scs=data_3d[index_scs_start-scs_shift:index_scs_end-scs_shift]*signal.tukey(len(scs_w_3d),alpha=0.5)
               scs_envp_3d = obsd_envp[index_scs_start-scs_shift:index_scs_end-scs_shift]*signal.tukey(len(scs_w_3d),alpha=0.5)
               index_1d_scs,scs_max_1d=max(enumerate(np.abs(scs_w_1d)),key=operator.itemgetter(1))
               index_3d_scs,scs_max_3d=max(enumerate(np.abs(shifted_3d_scs)),key=operator.itemgetter(1))
           # for i in range(len(sss_w_1d)):
               energy_scs=scipy.integrate.trapezoid((scs_w_1d-shifted_3d_scs)**2,dx=delta)
               scs1d_sqr=scipy.integrate.trapezoid(scs_w_1d*scs_w_1d,dx=delta)
               scs3d_sqr=scipy.integrate.trapezoid(shifted_3d_scs*shifted_3d_scs,dx=delta)
               scs1d_scs3d=scipy.integrate.trapezoid(scs_w_1d*shifted_3d_scs,dx=delta)
               scs_3d=scipy.integrate.trapezoid(shifted_3d_scs**2,dx=delta)
               norm_scs=scipy.integrate.trapezoid(scs_w_1d**2,dx=delta)
               A_scs=energy_scs/norm_scs
               R_scs=(1/len(scs_w_1d))*(scs_3d/norm_scs)
               amp_1d_scs=spec_amp(scs_w_1d,fr,delta)
               
               amp_3d_scs=spec_amp(shifted_3d_scs,fr,delta)               
               ### Envelope Misfit for SSS
               


               envp_mf_scs = np.log(max(scs_envp_3d)/max(scs_envp_1d))
               max_scs_envp_3d=np.argmax(scs_envp_3d)
               max_scs_envp_1d=np.argmax(scs_envp_1d)

               #envp_mf_scs = 0.5*scipy.integrate.trapezoid(lne**2)
               
               
           except AttributeError:
                  print("No ScS Phase detected")

           try:
               sdiff_window_start=st[1].stats.sac.t4-50.0
               sdiff_window_end=sdiff_window_start+100.0
               index_sdiff_start=int(sdiff_window_start/delta)
               index_sdiff_end=int(sdiff_window_end/delta)
               sdiff_w_1d=data_1d[index_sdiff_start:index_sdiff_end]
               sdiff_envp_1d=synt_envp[index_sdiff_start:index_sdiff_end]
               sdiff_w_1d=sdiff_w_1d*signal.tukey(len(sdiff_w_1d),alpha=0.5)            
               sdiff_w_3d=data_3d[index_sdiff_start:index_sdiff_end]*signal.tukey(len(sdiff_w_1d),alpha=0.5)
               sdiff_w_3d=sdiff_w_3d*signal.tukey(len(sdiff_w_3d),alpha=0.5)
               sdiff_corr=signal.correlate(sdiff_w_1d,sdiff_w_3d,method='auto',mode='full')
               max_ind=np.argmax(sdiff_corr)
               sdiff_shift=max_ind-len(sdiff_w_1d)+1
               shifted_3d_sdiff=data_3d[index_sdiff_start-sdiff_shift:index_sdiff_end-sdiff_shift]*signal.tukey(len(sdiff_w_3d),alpha=0.5)
               sdiff_envp_3d = obsd_envp[index_sdiff_start-sdiff_shift:index_sdiff_end-sdiff_shift]*signal.tukey(len(sdiff_w_3d),alpha=0.5)
               index_1d_sdiff,sdiff_max_1d=max(enumerate(np.abs(sdiff_w_1d)),key=operator.itemgetter(1))
               index_3d_sdiff,sdiff_max_3d=max(enumerate(np.abs(shifted_3d_sdiff)),key=operator.itemgetter(1))
           # for i in range(len(sss_w_1d)):
               energy_sdiff=scipy.integrate.trapezoid((sdiff_w_1d-shifted_3d_sdiff)**2,dx=delta)
               sdiff1d_sqr=scipy.integrate.trapezoid(sdiff_w_1d*sdiff_w_1d,dx=delta)
               sdiff3d_sqr=scipy.integrate.trapezoid(shifted_3d_sdiff*shifted_3d_sdiff,dx=delta)
               sdiff1d_sdiff3d=scipy.integrate.trapezoid(sdiff_w_1d*shifted_3d_sdiff,dx=delta)
               sdiff_3d=scipy.integrate.trapezoid(shifted_3d_sdiff**2,dx=delta)
               norm_sdiff=scipy.integrate.trapezoid(sdiff_w_1d**2,dx=delta)
               A_sdiff=energy_sdiff/norm_sdiff
               R_sdiff=(1/len(sdiff_w_1d))*(sdiff_3d/norm_sdiff)
               amp_1d_sdiff=spec_amp(sdiff_w_1d,fr,delta)
               
               amp_3d_sdiff=spec_amp(shifted_3d_sdiff,fr,delta)               
               ### Envelope Misfit for SSS
               

               envp_mf_sdiff = np.log(max(sdiff_envp_3d)/max(sdiff_envp_1d))
               max_sdiff_envp_3d=np.argmax(sdiff_envp_3d)
               max_sdiff_envp_1d=np.argmax(sdiff_envp_1d)

               #envp_mf_sdiff = 0.5*scipy.integrate.trapezoid(lne**2)
               
               
           except AttributeError:
                  print("No Sdiff Phase detected")
           ##Plotting the S, SS and SSS windows:
           fig,ax=plt.subplots(3,2,figsize=(20,12))
           if s_window_start!=0 and s_window_end!=0:
              ax[0,0].plot(time_1d[index_s_start:index_s_end],s_w_1d,label='1D isotropic prem',color='k')
              ax[0,0].plot(time_3d[index_s_start:index_s_end],shifted_3d_s,label='S40RTS',color='r')
              ax[0,0].plot(time_1d[index_s_start:index_s_end],s_envp_1d, label='PREM envelope',linestyle='--',color='k')
              ax[0,0].plot(time_3d[index_s_start:index_s_end],s_envp_3d,label='S40RTS envelope', linestyle='--',color='r')
              #ax[0,0].axvline(x=time_1d[index_s_start+index_1d_s],color='k',linestyle='-')
              #ax[0,0].axvline(x=s_arrival,color='b',linestyle='-')
              #ax[0,0].axvline(x=time_3d[index_s_start+index_3d_s],color='r',linestyle='--')
              ax[0,0].axvline(x=time_1d[index_s_start+max_s_envp_1d], color='k',linestyle='--')
              ax[0,0].axvline(x=time_3d[index_s_start+max_s_envp_3d],color='r', linestyle='--')
              ax[0,0].grid()
              ax[0,0].legend()
              #ax[0,0].title.set_text('S wave window(shift='+str(s_shift*delta)+' seconds)')
              #ax[0,1].plot(time_1d[min(index_s_start,index_s_start-s_shift):max(index_s_end,index_s_end-s_shift)],data_1d[min(index_s_start,index_s_start-s_shift):max(index_s_end,index_s_end-s_shift)],label='1D isotropic prem',color='k')
              #ax[0,1].plot(time_3d[min(index_s_start,index_s_start-s_shift):max(index_s_end,index_s_end-s_shift)],data_3d[min(index_s_start,index_s_start-s_shift):max(index_s_end,index_s_end-s_shift)],label='S40RTS',color='r')
              #ax[0,1].axvline(x=s_arrival,color='b',linestyle='-')
              #ax[0,1].grid()
              #ax[0,1].legend()
              ax[0,0].title.set_text('S wave window(shift='+str(s_shift*delta)+' seconds)')

           if ss_window_start!=0 and ss_window_end!=0:
              ax[1,0].plot(time_1d[index_ss_start:index_ss_end],ss_w_1d,label='1D isotropic prem',color='k')
              ax[1,0].plot(time_3d[index_ss_start:index_ss_end],shifted_3d_ss,label='S40RTS',color='r')
              ax[1,0].plot(time_1d[index_ss_start:index_ss_end],ss_envp_1d, label='PREM envelope',linestyle='--',color='k')
              ax[1,0].plot(time_3d[index_ss_start:index_ss_end],ss_envp_3d,label='S40RTS envelope', linestyle='--',color='r')
              #ax[1,0].axvline(x=time_1d[index_1d_ss+index_ss_start],color='k',linestyle='-')
              #ax[1,0].axvline(x=time_1d[,color='b',linestyle='-')
              #ax[1,0].axvline(x=time_3d[index_3d_ss+index_ss_start],color='r',linestyle='--')
              ax[1,0].axvline(x=time_1d[index_ss_start+max_ss_envp_1d], color='k',linestyle='--')
              ax[1,0].axvline(x=time_3d[index_ss_start+max_ss_envp_3d],color='r', linestyle='--')
              

              ax[1,0].grid()
              ax[1,0].legend()
              ax[1,0].title.set_text('SS phase window(shift='+str(ss_shift*delta)+' seconds)')       
              #ax[1,1].plot(time_1d[min(index_ss_start,index_ss_start-ss_shift):max(index_ss_end,index_ss_end-ss_shift)],data_1d[min(index_ss_start,index_ss_start-ss_shift):max(index_ss_end,index_ss_end-ss_shift)],label='1D isotropic prem',color='k')
              #ax[1,1].plot(time_3d[min(index_ss_start,index_ss_start-ss_shift):max(index_ss_end,index_ss_end-ss_shift)],data_3d[min(index_ss_start,index_ss_start-ss_shift):max(index_ss_end,index_ss_end-ss_shift)],label='S40RTS',color='r')
              #ax[1,1].axvline(x=ss_arrival,color='b',linestyle='-')
              #ax[1,1].grid()
              #ax[1,1].legend()
              #ax[1,1].title.set_text('SS wave window(shift='+str(ss_shift*delta)+' seconds)')
        
           if sss_window_start!=0 and sss_window_end!=0:
              ax[2,0].plot(time_1d[index_sss_start:index_sss_end],sss_w_1d,label='1D isotropic prem',color='k')
              ax[2,0].plot(time_3d[index_sss_start:index_sss_end],shifted_3d_sss,label='S40RTS',color='r')
              ax[2,0].plot(time_1d[index_sss_start:index_sss_end],sss_envp_1d, label='PREM envelope',linestyle='--',color='k')
              ax[2,0].plot(time_3d[index_sss_start:index_sss_end],sss_envp_3d,label='S40RTS envelope', linestyle='--',color='r')
              #ax[2,0].axvline(x=time_1d[index_1d_sss+index_sss_start],color='k',linestyle='-')
              #ax[2,0].axvline(x=sss_arrival,color='b',linestyle='-')
              #ax[2,0].axvline(x=time_3d[index_3d_sss+index_sss_start],color='r',linestyle='--')
              ax[2,0].axvline(x=time_1d[index_sss_start+max_sss_envp_1d], color='k',linestyle='--')
              ax[2,0].axvline(x=time_3d[index_sss_start+max_sss_envp_3d],color='r', linestyle='--')

              ax[2,0].grid()
              ax[2,0].legend()
              ax[2,0].title.set_text('SSS phase window(shift='+str(sss_shift*delta)+' seconds)')
              #ax[2,1].plot(time_1d[min(index_sss_start,index_sss_start-sss_shift):max(index_sss_end,index_sss_end-sss_shift)],data_1d[min(index_sss_start,index_sss_start-sss_shift):max(index_sss_end,index_sss_end-sss_shift)],label='1D isotropic prem',color='k')
              #ax[2,1].plot(time_3d[min(index_sss_start,index_sss_start-sss_shift):max(index_sss_end,index_sss_end-sss_shift)],data_3d[min(index_sss_start,index_sss_start-sss_shift):max(index_sss_end,index_sss_end-sss_shift)],label='S40RTS',color='r')
              #ax[2,1].axvline(x=sss_arrival,color='b',linestyle='-')
              #ax[2,1].grid()
              #ax[2,1].legend()
              #ax[2,1].title.set_text('SSS wave window(shift='+str(sss_shift*delta)+' seconds)')
           
           if scs_window_start != 0 and scs_window_end != 0:
              ax[0,1].plot(time_1d[index_scs_start:index_scs_end],scs_w_1d,label='1D isotropic prem',color='k')
              ax[0,1].plot(time_3d[index_scs_start:index_scs_end],shifted_3d_scs,label='S40RTS',color='r')
              ax[0,1].plot(time_1d[index_scs_start:index_scs_end],scs_envp_1d, label='PREM envelope',linestyle='--',color='k')
              ax[0,1].plot(time_3d[index_scs_start:index_scs_end],scs_envp_3d,label='S40RTS envelope', linestyle='--',color='r')
              #ax[3,0].axvline(x=time_1d[index_1d_scs+index_scs_start],color='k',linestyle='-')
              #ax[2,0].axvline(x=sss_arrival,color='b',linestyle='-')
              #ax[2,0].axvline(x=time_3d[index_3d_sss+index_sss_start],color='r',linestyle='--')
              ax[0,1].axvline(x=time_1d[index_scs_start+max_scs_envp_1d], color='k',linestyle='--')
              ax[0,1].axvline(x=time_3d[index_scs_start+max_scs_envp_3d],color='r', linestyle='--')


              ax[0,1].grid()
              ax[0,1].legend()
              ax[0,1].title.set_text('ScS phase window(shift='+str(scs_shift*delta)+' seconds)')
           if sdiff_window_start != 0 and sdiff_window_end != 0:
              ax[1,1].plot(time_1d[index_sdiff_start:index_sdiff_end],sdiff_w_1d,label='1D isotropic prem',color='k')
              ax[1,1].plot(time_3d[index_sdiff_start:index_sdiff_end],shifted_3d_sdiff,label='S40RTS',color='r')
              ax[1,1].plot(time_1d[index_sdiff_start:index_sdiff_end],sdiff_envp_1d, label='PREM envelope',linestyle='--',color='k')
              ax[1,1].plot(time_3d[index_sdiff_start:index_sdiff_end],sdiff_envp_3d,label='S40RTS envelope', linestyle='--',color='r')
              #ax[4,0].axvline(x=time_1d[index_1d_sdiff+index_sdiff_start],color='k',linestyle='-')
              #ax[2,0].axvline(x=sss_arrival,color='b',linestyle='-')
              #ax[2,0].axvline(x=time_3d[index_3d_sss+index_sss_start],color='r',linestyle='--')
              ax[1,1].axvline(x=time_1d[index_sdiff_start+max_sdiff_envp_1d], color='k',linestyle='--')
              ax[1,1].axvline(x=time_3d[index_sdiff_start+max_sdiff_envp_3d],color='r', linestyle='--')

              ax[1,1].grid()
              ax[1,1].legend()
              ax[1,1].title.set_text('Sdiff phase window(shift='+str(sdiff_shift*delta)+' seconds)')

              
           plt.savefig(dir_3d+'/../figures/'+str(stations[sta])+'_windows.png')
           plt.close() 
                #Plotting and calculating amplitudes for S phase
           ax=[0,0]
           fig = plt.figure(figsize=(15,10))
           gs = fig.add_gridspec(2, 1)
           ax[0]=fig.add_subplot(gs[0])
           ax[0].plot(time_1d,data_1d,label='1d_isotropic_prem',color='k',linewidth=0.5)
           ax[0].plot(time_3d,data_3d,label='3d_S40RTS',color='r',linewidth=0.5)
           try:
               ax[0].set_xlim(0,st[0].stats.sac.t2+200)
           except AttributeError:
                  ax[0].set_xlim(0,5400)
           try:
               ax[0].axvline(x=st[1].stats.sac.t0,color='k',label='S phase')

           except AttributeError:
                  print("No S phase detected, Epicentral Distance is "+str(st[0].stats.sac.gcarc))
       
        ##Plotting and Calculating amplitudes for SS phase

           try:
               ax[0].axvline(x=st[1].stats.sac.t1,color='b',label='SS phase')
           except AttributeError:
                  print("No SS phase detected,Epicentral Distance is "+str(st[0].stats.sac.gcarc))
            
        #Plotting amplitudes for SSS phase

           try:
               ax[0].axvline(x=st[1].stats.sac.t2,color='y',label='SSS phase')

           except AttributeError:
                  print("No SSS phase detected,Epicentral Distance is "+str(st[0].stats.sac.gcarc))
           try:
               ax[0].axvline(x=st[1].stats.sac.t3,color='r',label='ScS phase')

           except AttributeError:
                  print("No ScS phase detected,Epicentral Distance is "+str(st[0].stats.sac.gcarc))
           try:
               ax[0].axvline(x=st[1].stats.sac.t4,color='g',label='Sdiff phase')

           except AttributeError:
                  print("No Sdiff phase detected,Epicentral Distance is "+str(st[0].stats.sac.gcarc))

           ax[0].set_xlabel('Time')
           ax[0].legend()
           ax[1]=fig.add_subplot(gs[1],projection=ccrs.PlateCarree())
           ax[1].set_global()
           ax[1].coastlines()
           ax[1].add_feature(cfeature.LAND)

           ax[1].scatter(evlo, evla, 200, marker="*", color="red", transform=ccrs.PlateCarree(), zorder=100)
           ax[1].scatter(stlo, stla, 200, marker="v", color="blue", transform=ccrs.PlateCarree(), zorder=100)
           ax[1].plot([evlo, stlo],[evla,stla], "k", transform=ccrs.Geodetic())
        
           plt.savefig(dir_3d+'/../figures/'+str(stations[sta])+'.png')
           plt.close()

## Writing windows to new file
           #print(s_shift,s1d_sqr)
           if s_window_start!=0:    
              with open(dir_3d+'/../s_windows.txt','a') as txt:
                   txt.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(s_window_start)+'  '+str(s_window_end)+'   '+str(s_max_1d)+'  '+str(s_max_3d)+' '+str(s_shift*delta)+' '+str(A_s)+' '+str(R_s)+' '+str(s1d_sqr)+' '+str(s3d_sqr)+' '+str(s1d_s3d)+' '+str(baz)+' '+str(envp_mf_s)+'\n')
              with open(dir_3d+'/../s_amplitudes.txt','a') as wa:
                   wa.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(amp_1d_s)+' '+str(amp_3d_s)+' '+str(baz)+' '+str(envp_mf_s)+'\n')

           if ss_window_start!=0:
              with open(dir_3d+'/../ss_windows.txt','a') as wrt:
                 wrt.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(ss_window_start)+' '+str(ss_window_end)+'   '+str(ss_max_1d)+'  '+str(ss_max_3d)+' '+str(ss_shift*delta)+' '+str(A_ss)+' '+str(R_ss)+' '+str(ss1d_sqr)+' '+str(ss3d_sqr)+' '+str(ss1d_ss3d)+' '+str(baz)+' '+str(envp_mf_ss)+'\n')
              with open(dir_3d+'/../ss_amplitudes.txt','a') as wa:
                   wa.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(amp_1d_ss)+' '+str(amp_3d_ss)+' '+str(baz)+' '+str(envp_mf_ss)+'\n')

           if sss_window_start!=0:
              with open(dir_3d+'/../sss_windows.txt','a') as wrt:
                wrt.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(sss_window_start)+' '+str(sss_window_end)+'   '+str(sss_max_1d)+'  '+str(sss_max_3d)+' '+str(sss_shift*delta)+' '+str(A_sss)+' '+str(R_sss)+' '+str(sss1d_sqr)+' '+str(sss3d_sqr)+' '+str(sss1d_sss3d)  +' '+str(baz)+' '+str(envp_mf_sss)+'\n')
              with open(dir_3d+'/../sss_amplitudes.txt','a') as wa:
                   wa.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(amp_1d_sss)+' '+str(amp_3d_sss)+' '+str(baz)+' '+str(envp_mf_sss)+'\n')
           if scs_window_start!=0:
              with open(dir_3d+'/../scs_windows.txt','a') as wrt:
                wrt.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(scs_window_start)+' '+str(scs_window_end)+'   '+str(scs_max_1d)+'  '+str(scs_max_3d)+' '+str(scs_shift*delta)+' '+str(A_scs)+' '+str(R_scs)+' '+str(scs1d_sqr)+' '+str(scs3d_sqr)+' '+str(scs1d_scs3d)  +' '+str(baz)+' '+str(envp_mf_scs)+'\n')
              with open(dir_3d+'/../scs_amplitudes.txt','a') as wa:
                   wa.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(amp_1d_scs)+' '+str(amp_3d_scs)+' '+str(baz)+' '+str(envp_mf_scs)+'\n')
           if sdiff_window_start!=0:
              with open(dir_3d+'/../sdiff_windows.txt','a') as wrt:
                wrt.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(sdiff_window_start)+' '+str(sdiff_window_end)+'   '+str(sdiff_max_1d)+'  '+str(sdiff_max_3d)+' '+str(sdiff_shift*delta)+' '+str(A_sdiff)+' '+str(R_sdiff)+' '+str(sdiff1d_sqr)+' '+str(sdiff3d_sqr)+' '+str(sdiff1d_sdiff3d)  +' '+str(baz)+' '+str(envp_mf_sdiff)+'\n')
              with open(dir_3d+'/../sdiff_amplitudes.txt','a') as wa:
                   wa.write(str(stations[sta])+'  '+str(distance[sta])+'  '+str(evla)+' '+str(evlo)+' '+str(stla)+' '+str(stlo)+' '+str(evdp)+' '+str(amp_1d_sdiff)+' '+str(amp_3d_sdiff)+' '+str(baz)+' '+str(envp_mf_sdiff)+'\n')

              with open(dir_3d+'/../stations_done.txt','a') as txt2:
                   txt2.write(str(sta)+' '+str(stations[sta])+' '+str(distance[sta])+'\n')
        else:
             with open(dir_3d+'/../stations_out_of_range.txt','a') as wrt2:
                  wrt2.write(str(stations[sta])+'    '+str(distance[sta])+'\n') 
fi=[]
key='run00*'
for run in glob.glob(key):
    fi.append(run)

poolsize=5
pool=mp.Pool(poolsize)
pool.map(processing,fi)



