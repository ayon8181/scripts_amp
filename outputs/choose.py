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
import geopy.distance
from geographiclib.geodesic import Geodesic

labels=["PREM_3DC","S40RTS_3DC","S40RTS_1DC","GLAD_M25","Real_Data"]
file_names=["_prem_3D","_3D","_1D","_glad","_real"]

us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"]

ev_list=[]
with open("./../event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])

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


def plot(k,ph):
    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))]
    df_plot=df_plot[(df_plot["7"]>30)]
    for i,df in enumerate(file_names):
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+df])<12)]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+1)+df])>0.75)]
        df_plot=df_plot[df_plot[str(18+k*11+5)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+4)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+3)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+3)+df])<1.5)]
        df_plot=df_plot[df_plot[str(18+k*11+4)+df].notna()]
        df_plot=df_plot[(np.abs(np.log(df_plot[str(18+k*11+5)+df])<1.5))]
        df_plot=df_plot[df_plot[str(18+k*11+6)+df].notna()]
        df_plot=df_plot[(np.abs(df_plot[str(18+k*11+6)+df])<1.5)]
    df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    (df_plot["1"].str.split(".").str[0].isin(gsn))]
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
    print("Pacific Positive"+"\n")
    df_print = df_plot[(df_plot['mid_lon']<-120) & (df_plot['mid_lon']>-180) & (df_plot['mid_lat'] > -5) \
                       & (df_plot['mid_lat']<5) & (df_plot[str(18+k*11+6)+"_glad"]>0) & (df_plot[str(18+k*11+6)+"_3D"]>0)].iloc[0]
    
    print(str(df_print["0"])+ str(df_print["1"])+"\n")
    for i,df in enumerate(file_names):
        print(str((df_print[str(18+k*11+6)+df]))+"\n")
    print("Pacific Negative"+"\n")
    df_print = df_plot[(df_plot['mid_lon']<-120) & (df_plot['mid_lon']>-180) & (df_plot['mid_lat'] > -5) \
                       & (df_plot['mid_lat']<5) & (df_plot[str(18+k*11+6)+"_glad"]<0) & (df_plot[str(18+k*11+6)+"_3D"]<0)].iloc[0]
    
    print(str(df_print["0"])+ str(df_print["1"])+"\n")
    for i,df in enumerate(file_names):
        print(str((df_print[str(18+k*11+6)+df]))+"\n")
    print("North Atlantic"+"\n")
    df_print = df_plot[(df_plot['mid_lon']<55) & (df_plot['mid_lon']>45) & (df_plot['mid_lat'] > 0) \
                       & (df_plot['mid_lat']<70) & (df_plot[str(18+k*11+6)+"_glad"]<0) & (df_plot[str(18+k*11+6)+"_3D"]<0) & (df_plot[str(18+k*11+6)+"_1D"]>0)].iloc[0]
    
    print(str(df_print["0"])+ str(df_print["1"])+"\n")
    for i,df in enumerate(file_names):
        print(str((df_print[str(18+k*11+6)+df]))+"\n")

    print("Indian Ocean"+"\n")
    df_print = df_plot[(df_plot['mid_lon']>50) & (df_plot['mid_lon']<70) & (df_plot['mid_lat'] <-25) \
                       & (df_plot['mid_lat']>-35) & (df_plot[str(18+k*11+6)+"_glad"]<0) & (df_plot[str(18+k*11+6)+"_3D"]<0) &(df_plot[str(18+k*11+6)+"_1D"]>0)].iloc[0]
    
    print(str(df_print["0"])+ str(df_print["1"])+"\n")
    for i,df in enumerate(file_names):
        print(str((df_print[str(18+k*11+6)+df]))+"\n")
    print("African Plume"+"\n")
    df_print = df_plot[(df_plot['mid_lon']>0) & (df_plot['mid_lon']<45) & (df_plot['mid_lat'] <0) \
                       & (df_plot['mid_lat']>-5) & (df_plot[str(18+k*11+6)+"_glad"]<0) & (df_plot[str(18+k*11+6)+"_3D"]<0) & (df_plot[str(18+k*11+6)+"_1D"]>0)].iloc[0]
    
    print(str(df_print["0"])+ str(df_print["1"])+"\n")
    for i,df in enumerate(file_names):
        print(str((df_print[str(18+k*11+6)+df]))+"\n")
   
plot(4,"Sdiff")

    
    