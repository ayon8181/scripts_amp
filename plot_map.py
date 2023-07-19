import pygmt
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import multiprocessing as mp
from itertools import product
import os
def get_mid(evla,evlo,stla,stlo,dist):
    ##Give npts as odd numbers
    
    points=pygmt.project(center=[evlo,evla],endpoint=[stlo,stla],generate=2)
    npts=len(points)
    if len(points)%2==0:
       mp_ind=int((npts)/2)
       mp_ind_2=mp_ind-1
       lon=points.r.to_numpy() 
       mp_lon=(lon[mp_ind]+lon[mp_ind_2])/2.0
       lat=points.s.to_numpy()
       mp_lat=(lat[mp_ind]+lat[mp_ind_2])/2.0
       if len(points)%3==0:
          third=int(npts/3.0)
          two_third=2*int(npts/3.0)

          third_lon=(lon[third]*2+lon[third-1])/3.0
          two_third_lon=(2*lon[two_third-1]+lon[two_third])/3.0
          
          third_lat=(lat[third]*2+lat[third-1])/3.0
          two_third_lat=(2*lat[two_third-1]+lat[two_third])/3.0
          lons=np.zeros(len(lon[third:two_third])+2)
          lats=np.zeros(len(lon[third:two_third])+2)
          lons[0]=third_lon
          lons[1:-1]=lon[third:two_third]
          lons[-1]=two_third_lon
          lats[0]=third_lat
          lats[1:-1]=lat[third:two_third]
          lats[-1]=two_third_lat
       elif len(points)%3==1:
          third=int(npts/3.0)
          two_third=2*int(npts/3.0)

          third_lon=lon[third]
          two_third_lon=lon[two_third]
          
          third_lat=lat[third]
          two_third_lat=lat[two_third]
          lons=np.zeros(len(lon[third:two_third+1]))
          lats=np.zeros(len(lon[third:two_third+1]))
          #lons.append(third_lon)
          lons=lon[third:two_third+1]
          #lons.append(two_third_lon)
          #lats.append(third_lat)
          lats=lat[third:two_third+1]
         # lats.append(two_third_lat)
       else:
          third=int(npts/3.0)
          two_third=2*int(npts/3.0)

          third_lon=(lon[third]+2*lon[third+1])/3.0
          two_third_lon=(2*lon[two_third]+lon[two_third+1])/3.0
          
          third_lat=(lat[third]+2*lat[third+1])/3.0
          two_third_lat=(2*lat[two_third]+lat[two_third+1])/3.0
          lons=np.zeros(len(lon[third:two_third])+2)
          lats=np.zeros(len(lon[third:two_third])+2)
          lons[0]=third_lon
          lons[1:-1]=lon[third+1:two_third+1]
          lons[-1]=two_third_lon
          lats[0]=third_lat
          lats[1:-1]=lat[third+1:two_third+1]
          lats[-1]=two_third_lat
       
       return [[mp_lon,mp_lat],lons,lats]
    else:
         mp_ind=int((npts-1)/2)
         lon=points.r.to_numpy()
         mp_lon=lon[mp_ind]
         lat=points.s.to_numpy()
         mp_lat=lat[mp_ind]
         if len(points)%3==0:
            third=int(npts/3.0)
            two_third=2*int(npts/3.0)

            third_lon=(lon[third]*2+lon[third-1])/3.0
            two_third_lon=(2*lon[two_third-1]+lon[two_third])/3.0
          
            third_lat=(lat[third]*2+lat[third-1])/3.0
            two_third_lat=(2*lat[two_third-1]+lat[two_third])/3.0
            lons=np.zeros(len(lon[third:two_third])+2)
            lats=np.zeros(len(lon[third:two_third])+2)
            lons[0]=third_lon
            lons[1:-1]=lon[third:two_third]
            lons[-1]=two_third_lon
            lats[0]=third_lat
            lats[1:-1]=lat[third:two_third]
            lats[-1]=two_third_lat
         elif len(points)%3==1:
              third=int(npts/3.0)
              two_third=2*int(npts/3.0)

              third_lon=lon[third]
              two_third_lon=lon[two_third]
          
              third_lat=lat[third]
              two_third_lat=lat[two_third]
              lons=np.zeros(len(lon[third:two_third+1]))
              lats=np.zeros(len(lon[third:two_third+1]))
          #lons.append(third_lon)
              lons=lon[third:two_third+1]
          #lons.append(two_third_lon)
          #lats.append(third_lat)
              lats=lat[third:two_third+1]
         else:
              third=int(npts/3.0)
              two_third=2*int(npts/3.0)

              third_lon=(lon[third]+2*lon[third+1])/3.0
              two_third_lon=(2*lon[two_third]+lon[two_third+1])/3.0
          
              third_lat=(lat[third]+2*lat[third+1])/3.0
              two_third_lat=(2*lat[two_third]+lat[two_third+1])/3.0
              lons=np.zeros(len(lon[third:two_third])+2)
              lats=np.zeros(len(lon[third:two_third])+2)
              lons[0]=third_lon
              lons[1:-1]=lon[third+1:two_third+1]
              lons[-1]=two_third_lon
              lats[0]=third_lat
              lats[1:-1]=lat[third+1:two_third+1]
              lats[-1]=two_third_lat
       
         return [[mp_lon,mp_lat],lons,lats]             
         
dist=[]
evla=[]
evlo=[]
stla=[]
stlo=[]
sss_ss=[]
ss_s=[]
sss_s=[]
s_shift=[]
ss_shift=[]
sss_shift=[]
s_3d1d=[]
ss_3d1d=[]
sss_3d1d=[]
scs_3d1d=[]
sdiff_3d1d=[]
baz=[]
leng=[]
envp_sss=[]

plot_dir="./plot_shallow_1D"
os.chdir("/home/ayon/WORK/new_events/synthetics/")
#Reading Values
with open('./final_ratio_shallow_1D.txt','r') as txt:
     reader=csv.reader(txt,delimiter=' ')
     for row in reader:
         if float(row[1])>30.0:
            dist.append(float(row[1]))
            #leng.append(float(row[1])/2.0)
            evla.append(float(row[2]))
            evlo.append(float(row[3]))
            stla.append(float(row[4]))
            stlo.append(float(row[5]))
            s_shift.append(float(row[7]))
            ss_shift.append(float(row[8]))
            sss_shift.append(float(row[9]))
            ss_s.append(float(row[10]))
            sss_ss.append(float(row[12]))
            sss_s.append(float(row[14]))
            s_3d1d.append(float(row[16]))
            ss_3d1d.append(float(row[17]))
            sss_3d1d.append(float(row[18]))
            baz.append(float(row[19]))
            envp_sss.append(float(row[20]))
            scs_3d1d.append(float(row[23]))
            sdiff_3d1d.append(float(row[24]))
            leng.append(0.2)
midpoints_lat=[]
midpoints_lon=[]
third_lon=[]
third_lat=[]
ss_s=np.log(ss_s)
sss_ss=np.log(sss_ss)
sss_s=np.log(sss_s)


#print(stlo)
#with open('S40RTS_1000.txt','r') as txt:
#     reader=csv.reader(txt,skipinitialspace=True,delimiter=' ')
#     for row in reader:
#         x.append(float(row[0]))
#         y.append(float(row[1]))
 #        z.append(float(row[2])) 
     
#print(x[0])
data2=pd.DataFrame({'stla':stla,'stlo':stlo})
data3=pd.DataFrame({'evla':evla,'evlo':evlo})
#print(min(data2.stlo))

print('done')
n=len(evla)
#paramlist=list(product(evla,evlo,stla,stlo,dist))
print('done')

poolsize=4
print(poolsize)
pool=mp.Pool(poolsize)
results=pool.starmap(get_mid,((evla[i],evlo[i],stla[i],stlo[i],dist[i]) for i in range(n)))
pool.close()
pool.join()

#print(results)

for i in range(len(results)):  
    midpoints_lat.append(results[i][0][1])
    midpoints_lon.append(results[i][0][0])
    third_lat.append(results[i][2])
    third_lon.append(results[i][1])

#for k in range(len(evla)):
#    lon,lat=get_mid(evla[k],evlo[k],stla[k],stlo[k],dist[k])
#    midpoints_lat.append(lat)
#    midpoints_lon.append(lon)
    #print(lon)
#print(midpoints_lat,midpoints_lon)
### Start figure for s_shift

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'s_shift':s_shift,'third_lon':third_lon, 'third_lat':third_lat,'baz': baz,'leng':leng})
#data_l=data[data["s_shift"]<=0]
#data_g=data[data["s_shift"]>0]
#data=pd.DataFrame({'mid_lon':evlo,'mid_lat':evla,'s_shift':s_shift,'baz':baz,'leng':leng})
data_l=data[data["s_shift"] <= 0]
data_g=data[data["s_shift"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-5,5,5])
#fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.s_shift,cmap=True,transparency=40)
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#for j in range(len(stlo)):
#    fig.plot(x=stlo[j],y=stla[j],style="t0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/s_shift_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-5,5,5])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.s_shift,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/s_shift.png")
plt.close()

sns.histplot(data=data, x='s_shift',bins=19,binrange=(-5,5))
plt.savefig(plot_dir+"/s_shift_hist.png")

## SS Shift

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'ss_shift':ss_shift,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data=data.query("`ss_shift`>=-15 and `ss_shift`<=15")
data_l=data[data["ss_shift"] <= 0]
data_g=data[data["ss_shift"] > 0]
#print(data.third_lat)
#print(data.third_lat[3],data.third_lon[3])
#print(data.mid_lon[3],data.mid_lat[3])
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',projection="W180/12c",cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-20,20,20],continuous=True)

#for j in range(len(third_lon)):
  #  fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="0.1p,red")
    

fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/ss_shift_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-5,5,5])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.ss_shift,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/ss_shift.png")
plt.close()

sns.histplot(data=data, x='ss_shift',bins=19,binrange=(-15,15))
plt.savefig(plot_dir+'/ss_shift_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'sss_shift':sss_shift,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data=data.query("`sss_shift`>=-15 and `sss_shift`<=15")
data_l=data[data["sss_shift"] <= 0]
data_g=data[data["sss_shift"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-15,15,15],continuous=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/sss_shift_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-5,5,5])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.sss_shift,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/sss_shift.png")
plt.close()

sns.histplot(data=data, x='sss_shift',bins=19,binrange=(-15,15))
plt.savefig(plot_dir+'/sss_shift_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'ss_s':ss_s,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["ss_s"] <= 0]
data_g=data[data["ss_s"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],continuous=True, reverse=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],fill=data.ss_s,cmap=True)
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/ss_s_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.ss_s,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/ss_s.png")
plt.close()


sns.histplot(data=data, x='ss_s',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/ss_s_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'sss_ss':sss_ss,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["sss_ss"] <= 0]
data_g=data[data["sss_ss"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",frame=True,projection="W-150/12c")
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],reverse=True,continuous=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/sss_ss_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.sss_ss,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/sss_ss.png")
plt.close()

sns.histplot(data=data, x='sss_ss',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/sss_ss_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'sss_s':sss_s,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["sss_s"] <= 0]
data_g=data[data["sss_s"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],reverse=True,continuous=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/sss_s_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.sss_s,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/sss_s.png")
plt.close()

sns.histplot(data=data, x='sss_s',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/sss_s_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'s_3d1d':s_3d1d,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["s_3d1d"] <= 0]
data_g=data[data["s_3d1d"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],continuous=True,reverse=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/s_3d1d_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.s_3d1d,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/s_3d1d.png")
plt.close()

sns.histplot(data=data, x='s_3d1d',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/s_3d1d_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'ss_3d1d':ss_3d1d,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["ss_3d1d"] <= 0]
data_g=data[data["ss_3d1d"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],continuous=True, reverse=True)
#for j in range(len(third_lon)):
 #   fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/ss_3d1d_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.ss_3d1d,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/ss_3d1d.png")
plt.close()

sns.histplot(data=data, x='ss_3d1d',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/ss_3d1d_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'sss_3d1d':sss_3d1d,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["sss_3d1d"] <= 0]
data_g=data[data["sss_3d1d"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],continuous=True,reverse=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/sss_3d1d_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.sss_3d1d,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/sss_3d_1d.png")
plt.close()

sns.histplot(data=data, x='sss_3d1d',bins=9,binrange=(-1,1))
plt.savefig(plot_dir+'/sss_3d1d_hist.png')


pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'sss_envp3d1d':envp_sss,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["sss_envp3d1d"] <= 0]
data_g=data[data["sss_envp3d1d"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],continuous=True,reverse=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/sss_envp3d1d_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.sss_envp3d1d,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/sss_envp3d1d.png")
plt.close()


sns.histplot(data=data, x='sss_envp3d1d',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/sss_envp3d1d_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'scs_3d1d':scs_3d1d,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["scs_3d1d"] <= 0]
data_g=data[data["scs_3d1d"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],continuous=True,reverse=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/scs_3d1d_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.scs_3d1d,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/scs_3d1d.png")
plt.close()


sns.histplot(data=data, x='scs_3d1d',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/scs_3d1d_hist.png')

pygmt.makecpt(cmap='polar',reverse=True,series=[-2,2])
data=pd.DataFrame({'mid_lon':midpoints_lon,'mid_lat':midpoints_lat,'sdiff_3d1d':envp_sss,'third_lon':third_lon, 'third_lat':third_lat,'baz':baz,'leng':leng})
data_l=data[data["sdiff_3d1d"] <= 0]
data_g=data[data["sdiff_3d1d"] > 0]
fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c",frame=True)
#fig.grdimage(grid='S40RTS_600.nc',cmap=True)
fig.coast(shorelines=True, frame=True)
#pygmt.makecpt(cmap="polar",series=[-1,1,1],continuous=True,reverse=True)
#for j in range(len(third_lon)):
#    fig.plot(x=data.third_lon[j],y=data.third_lat[j],pen="1p,red")
fig.plot(x=data_l.mid_lon,y=data_l.mid_lat,style="v0.05c",direction=[data_l.baz,data_l.leng],fill="red",pen="red",transparency=70)
fig.plot(x=data_g.mid_lon,y=data_g.mid_lat,style="v0.05c",direction=[data_g.baz,data_g.leng],fill="blue",pen="blue",transparency=70)
#fig.plot(x=data3.evlo,y=data3.evla,style="a0.1c")
#fig.plot(x=data2.stlo,y=data2.stla,style="t0.1c")
#fig.colorbar()
fig.savefig(plot_dir+"/sdiff_3d1d_line.png")
plt.close()

fig=pygmt.Figure()
fig.basemap(region="d",projection="W-150/12c")
fig.coast(shorelines=True, frame=True)
pygmt.makecpt(cmap="polar",series=[-1,1,1])
fig.plot(x=data.mid_lon,y=data.mid_lat,style="c0.05c",fill=data.sdiff_3d1d,cmap=True,transparency=40)
fig.colorbar()
fig.savefig(plot_dir+"/sdiff_3d1d.png")
plt.close()


sns.histplot(data=data, x='sdiff_3d1d',bins=19,binrange=(-1,1))
plt.savefig(plot_dir+'/sdiff_3d1d_hist.png')



