#!/usr/bin/env python
import pyasdf
import obspy
import numpy as np
from scipy import integrate
from scipy.signal import hilbert
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy.signal.cross_correlation import correlate_template
from obspy.imaging import beachball
from obspy import Stream
from scipy import fft
from pyasdf import ASDFDataSet
from geographiclib.geodesic import Geodesic
import time
import sys
import math
import csv
if len(sys.argv) == 5:
    event       = str(sys.argv[1])
    obsd_tag    = str(sys.argv[2])
    obsd_fname  = str(sys.argv[3])
    f           = float(sys.argv[4])
else:
     print("Enter Input Information in the following order $event_name $asdf_tag $asdf_file $central_freq")

def rad_pat_sh(rake,dip,az,st,tk_of):
    y = (np.cos(rake)*np.cos(dip)*np.cos(tk_of)*np.sin(az-st) +
       np.cos(rake)*np.sin(dip)*np.sin(tk_of)*np.cos(2*(az-st)) +
       np.sin(rake)*np.cos(2*dip)*np.cos(tk_of)*np.cos(az-st) -
       0.5*np.sin(rake)*np.sin(2*dip)*np.sin(tk_of)*np.sin(2*(az-st)))
    return y

def nodal_plane_filter(event):
    mt    = event.focal_mechanisms[0].moment_tensor.tensor
    #nod_p = beachball.mt2plane(beachball.MomentTensor(0,1,-1,0,0,0,26))
    nod_p = beachball.mt2plane(beachball.MomentTensor(mt.m_rr,mt.m_tt,mt.m_pp,mt.m_rt,mt.m_rp,mt.m_tp,0))
    strike= nod_p.strike
    dip   = nod_p.dip
    rake  = nod_p.rake
    l=0
    np_0=[]
    az=[]
    np_1=[]
    for i in range(1,360):
        rp=np.abs(rad_pat_sh(np.deg2rad(rake),np.deg2rad(dip),np.deg2rad(i),np.deg2rad(strike),np.deg2rad(20.0)))
        if rp<=0.001:
           np_0.append(i)       
    return np_0

def plane(nodal):
    if nodal<0:
       y=360+nodal
    elif nodal>=360:
         y=0+(nodal-360)
    else:
         y=nodal
    
    return y

def SNR(noise,data):
    n_rate    = noise[1]-noise[0]
    d_rate    = data[1]-data[0]
    n_power   = integrate.trapezoid(noise*noise,dx=n_rate)/(len(noise)*n_rate)
    d_power   = integrate.trapezoid(data*data,dx=d_rate)/(len(data)*d_rate)
    rat_power = d_power/n_power

    m_noise   = max(np.abs(noise))
    m_data    = max(np.abs(data))

    rat_max   = m_data/m_noise

    return rat_power, rat_max

def make_measure(obsd,synt,stationXML,phase_list,event_loc,obsd_tag,nod_p,noise_acc,f):                 #event_loc=[evla,evlo,evdp]
    """
    Error Code:
    1: Rejected as the station is on a nodal plane
    2: Fails Global Signal to Noise Ratio Filters

    Error for each phase:
    0: Succesfully measured windows
    1: No arrival for that phase
    2: Signal to Noise Ratio Fails for that phase
    3: Correlation Value Error
    4: Time Shift greater than 15 seconds
    5: Signal to Noise Ratio Fails for time_shift windows
    """
    if f != 0:
       f_min = 1/f-(1/f)*0.15
       f_max = 1/f+(1/f)*0.15
    elif f == 0:
        f_min = 1/50.0
        f_max = 1/20.0
    

    model=TauPyModel(model="prem")
    #create output arrays
    max_cc       = np.full(len(phase_list),np.nan)
    amp_rat_3d1d = np.full(len(phase_list),np.nan)
    env_rat_3d1d = np.full(len(phase_list),np.nan)
    time_shifts  = np.full(len(phase_list),np.nan)
    obsd_obsd    = np.full(len(phase_list),np.nan)
    synt_synt    = np.full(len(phase_list),np.nan)
    obsd_synt    = np.full(len(phase_list),np.nan)
    A1           = np.full(len(phase_list),np.nan)
    A2           = np.full(len(phase_list),np.nan)
    arrival_times= np.full(len(phase_list),np.nan)
    success      = np.full(len(phase_list),np.nan)

    ## Read the station and event information
    sta    = stationXML
    stla   = sta[0][0].latitude
    stlo   = sta[0][0].longitude
    stel   = sta[0][0].elevation
    name   = sta[0].code+sta[0][0].code
    evla   = event_loc[0]
    evlo   = event_loc[1]
    evdp   = event_loc[2]

    

    gcarc  = locations2degrees(stla, stlo, evla, evlo)
    geod   = Geodesic.WGS84.Inverse(evla,evlo,stla,stlo)
    az     = geod["azi1"]

    ##Nodal Plane Filtering
    for p in nod_p:
        if plane(az) >= plane(p-7.5) and plane(az) <= plane(p+7.5):
            
            return [1,name]


    if gcarc >= 10:
    
        ## get the waveforms
        try:

            obsd_data     = obsd.select(component="T")[0]
            synt_data     = synt.select(component="T")[0]
            
            obsd_data     = obsd_data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
            synt_data     = synt_data.filter("bandpass",freqmin=f_min,freqmax=f_max,corners=4,zerophase=True)
            sampling_rate = synt_data.stats.sampling_rate
            dt            = synt_data.stats.delta
            obsd_anltc    = hilbert(obsd_data.data)
            synt_anltc    = hilbert(synt_data.data)
            obsd_envp     = np.abs(obsd_anltc)
            synt_envp     = np.abs(synt_anltc)

            
            
            arrival = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=gcarc,phase_list=["P"])
            n_start_index=np.nan
            if len(arrival) != 0:
               arrival_time = arrival[0].time

               arr_index    = int(arrival_time*sampling_rate)
               n_start_index= int(arr_index-125*sampling_rate)
               n_end_index  = int(arr_index-25*sampling_rate)
    
            elif len(arrival) == 0:
                 arrival = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=gcarc,phase_list=["Pdiff"])
                 if len(arrival) != 0:
                    arrival_time = arrival[0].time

                    arr_index    = int(arrival_time*sampling_rate)
                    n_start_index= int(arr_index-125*sampling_rate)
                    n_end_index  = int(arr_index-25*sampling_rate)
            if n_start_index!=n_start_index:
                 return [3,name]
            
            o_noise  = obsd_data.data[n_start_index:n_end_index]
            s_noise  = synt_data.data[n_start_index:n_end_index]

            ## Check Global Data Quality

            g_i, g_m = SNR(o_noise,obsd_data.data)

            if g_i <= noise_acc[0] and g_m <= noise_acc[1]:
               
               return [2,name]



            for i,phase in enumerate(phase_list):
                arrivals = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=gcarc,phase_list=[phase])
                if len(arrivals) != 0:
                    arrival_time = arrivals[0].time           #arrival_time in seconds
                
                    ## Cut 100s windows central to the 
                    arr_index    = int(arrival_time*sampling_rate)
                    start_index  = int(arr_index-50*sampling_rate)
                    end_index    = int(arr_index+50*sampling_rate)

                    synt_wdow    = synt_data.data[start_index:end_index]
                    en_synt_wdow = synt_envp[start_index:end_index]
                    obsd_cut     = obsd_data.data[int(start_index-50*sampling_rate):int(end_index+50*sampling_rate)]


                    ## Filter Data on SNR
                    w_i, w_m = SNR(o_noise,obsd_cut)
                    if w_i <= noise_acc[2] and w_m <= noise_acc[3]:
                        success[i] = 2
                        continue 

                    ## Calculate cross-correlation co-efficients
                    try:
                        cc           = correlate_template(obsd_cut,synt_wdow,mode="valid",normalize="full",method="auto")
                    except ValueError:
                        success[i] = 3
                        continue
                    max_ind      = np.argmax(cc)
                    max_cc_val   = cc[max_ind]
                    shift        = int(max_ind - int(len(cc)/2)-1)
                    time_shift   = shift/sampling_rate

                    ## Filter out based on cross-correlation value
                    if np.abs(time_shift) >= 15:
                       success[i] = 4
                       continue

                    obsd_wdow    = obsd_data[start_index+shift:end_index+shift]
                    en_obsd_wdow = obsd_envp[start_index+shift:end_index+shift]

                    w_i, w_m = SNR(o_noise,obsd_wdow)
                    if w_i <= noise_acc[2] and w_m <= noise_acc[3]:
                        success[i] = 5
                        continue 
                    

                    ## Make Measurements on the two windows
                    sqrd_obsd    = integrate.trapz(obsd_wdow*obsd_wdow, dx=dt)
                    sqrd_synt    = integrate.trapz(synt_wdow*synt_wdow, dx=dt)
                    obsd_synt_val= integrate.trapz(obsd_wdow*synt_wdow, dx=dt)
                    
                    energy_ratio = sqrd_obsd/sqrd_synt

                    max_envp_obs = max(en_obsd_wdow)
                    max_envp_syn = max(en_synt_wdow)
                    

                    envp_ratio   = max_envp_obs/max_envp_syn

                    max_cc[i]       = max_cc_val
                    amp_rat_3d1d[i] = energy_ratio
                    env_rat_3d1d[i] = envp_ratio
                    time_shifts[i]  = time_shift
                    obsd_obsd[i]    = sqrd_obsd
                    synt_synt[i]    = sqrd_synt
                    obsd_synt[i]    = obsd_synt_val
                    arrival_times[i]= arrival_time

                    A1[i]           = obsd_synt_val/sqrd_synt
                    A2[i]           = sqrd_obsd/obsd_synt_val
                    success[i]      = 0
                else:
                     success[i] = 1
            
            SS_S            = 0.5*(min(A1[1],A2[1])/max(A1[0],A2[0]) + max(A1[1],A2[1])/min(A1[0],A2[0]))
            SSS_SS          = 0.5*(min(A1[2],A2[2])/max(A1[1],A2[1]) + max(A1[2],A2[2])/min(A1[1],A2[1]))
            SSS_ScS         = 0.5*(min(A1[2],A2[2])/max(A1[3],A2[3]) + max(A1[2],A2[2])/min(A1[3],A2[3]))
            SSS_Sdiff       = 0.5*(min(A1[2],A2[2])/max(A1[4],A2[4]) + max(A1[2],A2[2])/min(A1[4],A2[4]))

            SS_S_error      = 0.5*(min(A1[1],A2[1])/max(A1[0],A2[0]) - max(A1[1],A2[1])/min(A1[0],A2[0]))
            SSS_SS_error    = 0.5*(min(A1[2],A2[2])/max(A1[1],A2[1]) - max(A1[2],A2[2])/min(A1[1],A2[1]))
            SSS_ScS_error   = 0.5*(min(A1[2],A2[2])/max(A1[3],A2[3]) - max(A1[2],A2[2])/min(A1[3],A2[3]))
            SSS_Sdiff_error = 0.5*(min(A1[2],A2[2])/max(A1[4],A2[4]) - max(A1[2],A2[2])/min(A1[4],A2[4]))

            output=[evlo, evla, evdp, stlo, stla, gcarc, SS_S, SS_S_error, SSS_SS, SSS_SS_error, SSS_ScS, SSS_ScS_error, SSS_Sdiff, SSS_Sdiff_error]
            for j,k in enumerate(phase_list):
                output.append(arrival_times[j])
                output.append(max_cc[j])
                output.append(time_shifts[j])
                output.append(amp_rat_3d1d[j])
                output.append(env_rat_3d1d[j])
                output.append(0.5*(A1[j]+A2[j]))
            
            return [0,name,output,success]
        except IndexError:
               return [3,name]





ds        = ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s.proc_"+obsd_fname+".h5",mode="r")
ds2       = ASDFDataSet("/scratch1/09038/ayon8181/scripts_amp/seis/proc/"+event+".T017-250s.proc_synt.h5",mode="r")

noise_acc = [3.5,6.0,3,3]
#print("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s.proc_"+obsd_fname+".h5")

evla      = ds.events[0].origins[0].latitude
evlo      = ds.events[0].origins[0].longitude
evdp      = ds.events[0].origins[0].depth/1000.0
evloc     = [evla,evlo,evdp]
phase_list= ["S","SS","SSS","ScS","Sdiff"] 

nod_p     = nodal_plane_filter(ds.events[0])


def process(this_station_group, other_station_group):
        # Make sure everything thats required is there.
    #if (not hasattr(this_station_group, "StationXML")
    #    or not hasattr(this_station_group, "proc_"+obsd_tag)
    #    or not hasattr(other_station_group, "proc_synt")):
    
    #    return
    

    stationxml = this_station_group.StationXML
    observed   = this_station_group["proc_"+obsd_tag]
    synthetic  = other_station_group.proc_synt

    all_results= []

    results    = make_measure(observed,synthetic,stationxml,phase_list,evloc,obsd_tag,nod_p,noise_acc,f)


    all_results.append(results)

    return all_results


a          = time.time
all_output = ds.process_two_files(ds2, process)
b          = time.time



if ds.mpi.rank == 0:
    #print(all_output)
    with open("/scratch1/09038/ayon8181/scripts_amp/outputs/"+str(f)+"_"+obsd_fname+".txt",'a') as txt, \
         open("/scratch1/09038/ayon8181/scripts_amp/errors/error_"+str(f)+"_"+obsd_fname+"_phase.txt",'a') as txt2, \
         open("/scratch1/09038/ayon8181/scripts_amp/errors/error_"+str(f)+"_"+obsd_fname+"_all.txt",'a') as txt3:
        for k in all_output.keys():
            if all_output[k][0] is not None and all_output[k][0][0] == 0 :
                txt.write(event+" "+k+" ")
                for h in all_output[k][0][2]:
                    txt.write(str(h)+" ")
                txt.write("\n")

                for j,error in enumerate(all_output[k][0][3]):
                    if error is not np.nan and error != 0:
                        txt2.write(event+" "+str(all_output[k][0][1])+" "+phase_list[j]+" "+str(error)+"\n")
                

            elif all_output[k][0] is not None and all_output[k][0][0] != 0:
                
                txt3.write(event+" "+str(all_output[k][0][1])+" "+str(all_output[k][0][0])+"/n")
            
del ds
del ds2    







    


        









