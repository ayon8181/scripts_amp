import pygmt
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import product
import os
import glob
import obspy

home="/home1/09038/ayon8181/"
lats=[]
lons=[]
depth=[]
mags=[]
gsn=["II","IU","CU","IC","GT","US","CN"]
us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
stations={}
events=[]
ev  =[home+"deep_100",home+"shallow_100"]
for i in ev:
    os.chdir(i)
    for j in glob.glob("*.txt"):
        st=obspy.read_events(j)
        evdp=st[0].origins[0].depth
        evla=st[0].origins[0].latitude
        evlo=st[0].origins[0].longitude
        mag =st[0].magnitudes[0].mag
        #events.append(j.split(".")[0])
        depth.append(evdp/1000.0)
        lats.append(evla)
        lons.append(evlo)
        mags.append(0.05*mag)
    
    for k in glob.glob("STATION_*"):
        with open(k,"r") as txt:
             reader=csv.reader(txt,skipinitialspace=True,delimiter="\t")
             for row in reader:
                
                name=row[1]+"."+row[0]
                stla=float(row[2])
                stlo=float(row[3])
                #if ((((stlo>=-135 and stlo<=-59) and (stla>=32 and stla<=80)) or ((stlo>=-166 and stlo<=-132) and (stla>=54 and stla<=71)) or (row[1] in us) or ((stlo>=-103 and stlo<=-59) and (stla>=25 and stla<=32))) and (row[1] not in gsn)): 
                   #f=1               
                   
                if name not in stations.keys():
                    stations[name]=[stla,stlo]

stlas=[]
stlos=[]
for l in stations.keys():
    stlas.append(stations[l][0])
    stlos.append(stations[l][1])
#rint(stlas)
#print(stlos)
#
max_depth=max(depth)
min_depth=min(depth)
fig=pygmt.Figure()
fig.basemap(region="d",projection="N-150/12c",frame=True)

las=[]
los=[]
with open("/scratch1/09038/ayon8181/scripts_amp/outputs/STATIONS","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="\t")
     
     for row in reader:
        # print(row)
         las.append(float(row[2]))
         los.append(float(row[4]))
# define etopo data file

fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="jet",series=[0,650,50])

fig.plot(x=stlos,y=stlas,style = 't0.05c',color = 'black',pen = 'thinner,black')
fig.plot(x=los,y=las,style = 't0.05c',color = 'red',pen = 'thinner,red')
fig.plot(x=lons,y=lats,size=mags,fill=depth,cmap=True,style="cc",pen="black")
fig.colorbar(frame=["xa100f0.1+lDepth", "y+lkm"])
fig.savefig("/scratch1/09038/ayon8181/scripts_amp/sou_sta.png")

        
        
