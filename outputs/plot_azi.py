import pandas as pd
import matplotlibpyplot as plt
import csv
from geographiclib.geodesic import Geodesic


g=Geodesic.WGS84
evla=[]
evlo=[]
stla=[]
stlo=[]
value=[]
with open("/scratch1/09038/ayon8181/synthetic/scripts_amp/outputs/figures/event/201906281551A/amp_file__glad_Sdiff.csv") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter=" ")
     next(reader)
     for row in reader:
         evla.append(float(row[3]))
         evlo.append(float(row[2]))
         stla.append(float(row[6]))
         stlo.append(float(row[5]))
         value.append(float(row[-1]))
az=[]
for i,k in enumerate(evla):
    az.append(g.Inverse(evla[i],evlo[i],stla[i],stlo[i]['azi1']))

plt.plot(az,value)
plt.show()



