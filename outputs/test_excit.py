import csv
import math
import obspy
import glob
from obspy.core.event import read_events
moho        = {}
with open("/scratch1/09038/ayon8181/scripts_amp/outputs/depthtomoho.xyz",'r') as txt:
    data = csv.reader(txt, skipinitialspace=True, delimiter=" ")
    for row in data:
        lats   = math.floor(float(row[1]))
        lons   = math.floor(float(row[0]))
        if lats not in moho.keys():
            moho[lats] = {}
        moho[lats][lons] = float(row[2])

for i in glob.glob("/work2/09038/ayon8181/frontera/new_events/synthetics/quakeml/*xml"):
    cat=read_events(i)
    print(cat[0].event_descriptions)
    depth=cat[0].origins[0].depth
    lat=math.floor(cat[0].origins[0].latitude)
    lon=math.floor(cat[0].origins[0].longitude)
    i_1=0
    i_2=0
    if -depth/1000 > moho[lat][lon]:
        print("Crust in 3D Crust")
        i_1=1 
    else:
         print("Mantle in 3D Crust") 
         i_1=0
    if -depth/1000 > -24.4:
        print("Crust in 1D Crust") 
        i_2=1
    else:
         print("Mantle in 1D Crust")
         i_2=0 
    if i_1 != i_2:
       with open("change.txt","a") as txt:

            #txt.write(str(cat[0].resource_id).split("/")[2]+" "+str(moho[lat][lon])+" "+str(depth)+"\n")
            txt.write(str(cat[0].resource_id).split("/")[2]+" "+str(depth)+" "+str(moho[lat][lon])+"\n")
    #else:
    #    with open("change.txt","a") as txt:

            #txt.write(str(cat[0].resource_id).split("/")[2]+" "+str(moho[lat][lon])+" "+str(depth)+"\n")
    #        txt.write(str(cat[0].resource_id).split("/")[2]+"\n")
print(moho[-70][-30])
