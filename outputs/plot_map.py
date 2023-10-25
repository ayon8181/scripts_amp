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

file_names      = ["obsd_1D_crust", "obsd_3D_crust","obsd_glad"]
df=[]
plot_dir="/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/plots_map"
model=TauPyModel(model="prem")
for i,f in enumerate(file_names):
    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt4",header=None,delimiter=" ",skipinitialspace=True)
    df.append(temp)

def points_S_ScS(evla, evlo, evdp, stla, stlo, model, phase_list):

    for j,p in enumerate(phase_list):
        k = model.get_ray_paths_geo(source_depth_in_km=evdp, source_latitude_in_deg=evla, source_longitude_in_deg=evlo, receiver_latitude_in_deg=stla, receiver_longitude_in_deg=stlo,phase_list=[p],resample=True) 
        if len(k) != 0:
           path=k[0].path
           mid=[]
           deep=0
           start=[path[int(len(path)/3)][5],path[int(len(path)/3)][4]]
           for i,pts in enumerate(path):
               if pts[3]>deep:
                  deep=pts[3]
                  mid=[pts[5],pts[4]]
           end=[path[int(len(path)*2/3)][5],path[int(len(path)*2/3)][4]]
           return [start,mid,end]

def points_SS(evla,evlo,evdp,stla,stlo,model=model):
    k = model.get_ray_paths_geo(source_depth_in_km=evdp, source_latitude_in_deg=evla, source_longitude_in_deg=evlo, receiver_latitude_in_deg=stla, receiver_longitude_in_deg=stlo,phase_list=['SS'],resample=True)
    path=k[0].path
    start=[]
    g=0
    if len(k) != 0:
        for i,pts in enumerate(path[0:-2]):
            if path[i][3]>path[i+1][3]:
                if len(start) !=2:
                    start=[pts[5],pts[4]]
            if pts[3] == 0.0:
                mid=[pts[5],pts[4]]
                g=1
                end=[]
            if path[i][3]>path[i+1][3] and g==1:
                if len(end) != 2:
                    end = [pts[5],pts[4]]
        return [start,mid,end]
    
def points_SSS(evla,evlo,evdp,stla,stlo,model=model):

    k = model.get_ray_paths_geo(source_depth_in_km=evdp, source_latitude_in_deg=evla, source_longitude_in_deg=evlo, receiver_latitude_in_deg=stla, receiver_longitude_in_deg=stlo,phase_list=['SSS'],resample=True)
    path=k[0].path
    start=[]
    if len(k) != 0:
        start=[]
        end=[]
        y=0
        deep=0
        for i,pts in enumerate(path[0:-2]):
            if pts[3] == 0.0:
                if len(start)==0:
                    start=[pts[5],pts[4]]
                    y=1
                elif len(end)==0:
                    end=[pts[5],pts[4]]
                    y=0
                if y==1:
                    if pts[3] >= deep:
                        mid=[pts[5],pts[4]]
                        deep = pts[3]
                     
                        
                  
                  
        return [start,mid,end]
def points_Sdiff(evla,evlo,evdp,stla,stlo,model=model):
    k = model.get_ray_paths_geo(source_depth_in_km=evdp, source_latitude_in_deg=evla, source_longitude_in_deg=evlo, receiver_latitude_in_deg=stla, receiver_longitude_in_deg=stlo,phase_list=['Sdiff'],resample=True) 
    if len(k) != 0:
        path=k[0].path
        mid=[]
        deep=0
        start=[path[int(len(path)/3)][5],path[int(len(path)/3)][4]]
        for i,pts in enumerate(path):
            if pts[3]==2190:
                mid.append([pts[5],pts[4]])
        mid=mid[int(len(mid)/2)] 
        end=[path[int(len(path)*2/3)][5],path[int(len(path)*2/3)][4]]
        return [start,mid,end]
                     
                        
                  
                  

phase_list=['S','SS','SSS','ScS','Sdiff']
for k,ph in enumerate(phase_list):
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[23+k*5+1]>0.85) & (np.abs((dfs[23+k*5+2])<15))]
        df_plot = df_plot[df_plot[23+k*5+3] == df_plot[23+k*5+3]]
        if k==0:
           df_plot = df_plot[(df_plot[6]>30) & (df_plot[6]<70)]#['0','1','2','3','4',str(23+k*5+3),str(23+k*5+4),str(23+k*5+2)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_S_ScS, [(rows[2],rows[1],rows[3],rows[5],rows[4],model,['S'])  for i,rows in df_plot.iterrows()])
           pool.close()
           pool.join()
           start=[]
           mid_lat=[]
           mid_lon=[]
           end=[]
           for r in results:
               if r is not None:
                  start.append(r[0]+r[1])
                  mid_lat.append(r[1][1])
                  mid_lon.append(r[1][0])
                  end.append(r[1]+r[2])
                  
               else:
                    start.append(np.nan)
                    
                    mid_lat.append(np.nan)
                    mid_lon.append(np.nan)
                    end.append(np.nan)
                    
           df_plot['start'] = start
           df_plot['mid_lon']   = mid_lon
           df_plot['mid_lat']   = mid_lat
           df_plot['end']   = end
        elif k==1:
           df_plot = df_plot[(df_plot[6]>50) & (df_plot[6]<140)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_SS, [(rows[2],rows[1],rows[3],rows[5],rows[4],model)  for i,rows in df_plot.iterrows()])
           pool.close()
           pool.join()
           start=[]
           mid_lat=[]
           mid_lon=[]
           end=[]
           for r in results:
               if r is not None:
                  start.append(r[0]+r[1])
                  mid_lat.append(r[1][1])
                  mid_lon.append(r[1][0])
                  end.append(r[1]+r[2])
                  
               else:
                    start.append(np.nan)
                    
                    mid_lat.append(np.nan)
                    mid_lon.append(np.nan)
                    end.append(np.nan)
                    
           df_plot['start'] = start
           df_plot['mid_lon']   = mid_lon
           df_plot['mid_lat']   = mid_lat
           df_plot['end']   = end
        elif k==2:
           df_plot = df_plot[((df_plot[6]>75) & (df_plot[6]<110)) | ((df_plot[6]>110) & (df_plot[6]<165))]#['0','1','2','3','4',str(23+k*5+3),str(23+k*5+4),str(23+k*5+2)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_SSS, [(rows[2],rows[1],rows[3],rows[5],rows[4],model)  for i,rows in df_plot.iterrows()])
           pool.close()
           pool.join()
           start=[]
           mid_lat=[]
           mid_lon=[]
           end=[]
           for r in results:
               if r is not None:
                  start.append(r[0]+r[1])
                  mid_lat.append(r[1][1])
                  mid_lon.append(r[1][0])
                  end.append(r[1]+r[2])
                  
               else:
                    start.append(np.nan)
                    
                    mid_lat.append(np.nan)
                    mid_lon.append(np.nan)
                    end.append(np.nan)
                    
           df_plot['start'] = start
           df_plot['mid_lon']   = mid_lon
           df_plot['mid_lat']   = mid_lat
           df_plot['end']   = end
        elif k==3:
           df_plot = df_plot[((df_plot[6]>5) & (df_plot[6]<25)) | ((df_plot[6]>50) & (df_plot[6]<65)) ]#['0','1','2','3','4',str(23+k*5+3),str(23+k*5+4),str(23+k*5+2)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_S_ScS, [(rows[2],rows[1],rows[3],rows[5],rows[4],model,['ScS'])  for i,rows in df_plot.iterrows()])
           pool.close()
           pool.join()
           start=[]
           mid_lat=[]
           mid_lon=[]
           end=[]
           for r in results:
               if r is not None:
                  start.append(r[0]+r[1])
                  mid_lat.append(r[1][1])
                  mid_lon.append(r[1][0])
                  end.append(r[1]+r[2])
                  
               else:
                    start.append(np.nan)
                    
                    mid_lat.append(np.nan)
                    mid_lon.append(np.nan)
                    end.append(np.nan)
                    
           df_plot['start'] = start
           df_plot['mid_lon']   = mid_lon
           df_plot['mid_lat']   = mid_lat
           df_plot['end']   = end
        elif k==4:
           df_plot = df_plot[(df_plot[6]>100) & (df_plot[6]<150)]['0','1','2','3','4',str(23+k*5+3),str(23+k*5+4),str(23+k*5+2)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_Sdiff, [(rows[2],rows[1],rows[3],rows[5],rows[4],model)  for i,rows in df_plot.iterrows()])
           pool.close()
           pool.join()
           start=[]
           mid_lat=[]
           mid_lon=[]
           end=[]
           for r in results:
               if r is not None:
                  start.append(r[0]+r[1])
                  mid_lat.append(r[1][1])
                  mid_lon.append(r[1][0])
                  end.append(r[1]+r[2])
                  
               else:
                    start.append(np.nan)
                    
                    mid_lat.append(np.nan)
                    mid_lon.append(np.nan)
                    end.append(np.nan)
                    
           df_plot['start'] = start
           df_plot['mid_lon']   = mid_lon
           df_plot['mid_lat']   = mid_lat
           df_plot['end']   = end
        print(df_plot)
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
        fig.coast(shorelines=True, frame=True)   
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.05c",fill=np.log(df_plot[23+k*5+3]),cmap=True,transparency=40)
        #fig.plot(data=df_plot.start,style="=0.0c",pen="1p,+z",zvalue=df_plot[23+k*5+3],cmap=True,transparency=40)            
        #fig.plot(data=df_plot.end,style="=0.0c",pen="1p,+z",zvalue=df_plot[23+k*5+3],cmap=True,transparency=40)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_amp"+file_names[i]+"_"+ph+".png")
        plt.close()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
        fig.coast(shorelines=True, frame=True)   
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.05c",fill=np.log(df_plot[23+k*5+4]),cmap=True,transparency=40)
        #fig.plot(data=df_plot.start,style="=0.0c",pen="1p,+z",zvalue=df_plot[23+k*5+4],cmap=True,transparency=40)            
        #fig.plot(data=df_plot.end,style="=0.0c",pen="1p,+z",zvalue=df_plot[23+k*5+4],cmap=True,transparency=40)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_env"+file_names[i]+"_"+ph+".png")
        plt.close()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-10,10]) 
        fig.basemap(region="d",projection="N-150/12c")
        fig.coast(shorelines=True, frame=True)   
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.05c",fill=df_plot[23+k*5+2],cmap=True,transparency=40)
        #fig.plot(data=df_plot.start,style="=0.0c",pen="1p,+z",zvalue=df_plot[23+k*5+2],cmap=True,transparency=40)            
        #fig.plot(data=df_plot.end,style="=0.0c",pen="1p,+z",zvalue=df_plot[23+k*5+2],cmap=True,transparency=40)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_time"+file_names[i]+"_"+ph+".png")
        plt.close
