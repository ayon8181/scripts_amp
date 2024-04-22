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

stations={}

ev  =["shallow_100","deep_100"]
for i in ev:
    os.chdir(home+i)
    for j in glob.glob("*.txt"):
        st=obspy.read_events(j)
        evdp=st[0].origins[0].depth
        evla=st[0].origins[0].latitude
        evlo=st[0].origins[0].longitude
        mag =st[0].magnitudes[0].mag

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
                 if name not in stations.keys():
                    stations[name]=[stla,stlo]

stlas=[]
stlos=[]
for l in stations.keys():
    stlas.append(stations[l][0])
    stlos.append(stations[l][1])

max_depth=max(depth)
min_depth=min(depth)
fig=pygmt.Figure()
fig.basemap(region="d",projection="N-150/12c",frame=True)


# define etopo data file

fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="jet",series=[0,650,50])

fig.plot(x=stlos,y=stlas,style = 't0.15c',color = 'black',pen = 'thinner,white')
fig.plot(x=lons,y=lats,size=mags,fill=depth,cmap=True,style="cc",pen="black")
fig.colorbar()
fig.savefig("/scratch1/09038/ayon8181/scripts_amp/outputs/sou_sta.png")

        
        
