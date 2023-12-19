#!/usr/bin/env python
import pyasdf
import obspy
import numpy as np
from scipy import integrate
from scipy.signal import hilbert
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy.signal.cross_correlation import correlate_template
from scipy import fft
from pyasdf import ASDFDataSet
import time
import sys
import math
import csv

event       = str(sys.argv[1])
obsd_tag    = str(sys.argv[2])
obsd_fname  = str(sys.argv[3])


def make_measure(obsd,synt,stationXML,phase_list,event_loc,obsd_tag):                 #event_loc=[evla,evlo,evdp]
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

    ## Read the station and event information
    sta    = stationXML
    stla   = sta[0][0].latitude
    stlo   = sta[0][0].longitude
    stel   = sta[0][0].elevation
    evla   = event_loc[0]
    evlo   = event_loc[1]
    evdp   = event_loc[2]
    

    gcarc= locations2degrees(stla, stlo, evla, evlo)
    if gcarc >= 10:
    
        ## get the waveforms
        try:
            obsd_data     = obsd.select(component="T")[0]
            synt_data     = synt.select(component="T")[0]
            sampling_rate = synt_data.stats.sampling_rate
            dt            = synt_data.stats.delta
            obsd_anltc    = hilbert(obsd_data.data)
            synt_anltc    = hilbert(synt_data.data)
            obsd_envp     = np.abs(obsd_anltc)
            synt_envp     = np.abs(synt_anltc)

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
                    try:
                        cc           = correlate_template(obsd_cut,synt_wdow,mode="valid",normalize="full",method="auto")
                    except ValueError:
                        return
                    max_ind      = np.argmax(cc)
                    max_cc_val   = cc[max_ind]
                    shift        = int(max_ind - int(len(cc)/2)-1)
                    time_shift   = shift/sampling_rate
                    obsd_wdow    = obsd_data[start_index+shift:end_index+shift]
                    en_obsd_wdow = obsd_envp[start_index+shift:end_index+shift]

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
            
            return output
        except IndexError:
               return None





ds        = ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-040s.proc_"+obsd_fname+".h5",mode="r")
ds2       = ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-040s.proc_real_data.h5",mode="r")
#print("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/"+event+".T017-050s.proc_"+obsd_fname+".h5")

evla      = ds.events[0].origins[0].latitude
evlo      = ds.events[0].origins[0].longitude
evdp      = ds.events[0].origins[0].depth/1000.0
evloc     = [evla,evlo,evdp]
phase_list= ["S","SS","SSS","ScS","Sdiff"] 


def process(this_station_group, other_station_group):
        # Make sure everything thats required is there.
    #if (not hasattr(this_station_group, "StationXML")
    #    or not hasattr(this_station_group, "proc_"+obsd_tag)
    #    or not hasattr(other_station_group, "proc_synt")):
    
    #    return
    

    stationxml = this_station_group.StationXML
    observed   = this_station_group.proc_real_data
    synthetic  = other_station_group["proc_"+obsd_tag]

    all_results= []

    results    = make_measure(observed,synthetic,stationxml,phase_list,evloc,obsd_tag)


    all_results.append(results)

    return all_results


a          = time.time
all_output = ds2.process_two_files(ds, process)
b          = time.time



if ds.mpi.rank == 0:
    #print(all_output)
    with open("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+obsd_fname+"_w_real.txt4",'a') as txt:
        for k in all_output.keys():
            if all_output[k][0] is not None:
                txt.write(event+" "+k+" ")
                for h in all_output[k][0]:
                    txt.write(str(h)+" ")
                txt.write("\n")
del ds
del ds2    







    


        









