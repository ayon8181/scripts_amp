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

file_names      = ["obsd_1D_crust", "obsd_3D_crust","obsd_glad","real_data"]
df=[]
plot_dir="/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/plots_map_deep"
model=TauPyModel(model="prem")
for i,f in enumerate(file_names):
    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt2",header=None,delimiter=" ",skipinitialspace=True)
    #temp=temp[temp[4]>=100]
    df.append(temp)
max_cc        = [0.0,0.0,0.0,0.0]
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
        if len(mid)>0:
           mid=mid[int(len(mid)/2)] 
           end=[path[int(len(path)*2/3)][5],path[int(len(path)*2/3)][4]]
           return [start,mid,end]
                     
                        
                  
                  

phase_list=['S','SS','SSS','ScS','Sdiff']
for k,ph in enumerate(phase_list):
    for i,dfs in enumerate(df):
        df_plot = dfs[(dfs[16+k*6+1]>max_cc[i]) & (np.abs((dfs[16+k*6+2])<15))]
        df_plot = df_plot[df_plot[16+k*6+3] == df_plot[16+k*6+3]]
        if k==0:
           df_plot = df_plot[(df_plot[7]>30) & (df_plot[7]<70)]#['0','1','2','3','4',str(16+k*6+3),str(16+k*6+4),str(16+k*6+2)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_S_ScS, [(rows[3],rows[2],rows[4],rows[6],rows[5],model,['S'])  for i,rows in df_plot.iterrows()])
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
        elif k==1:
           df_plot = df_plot[(df_plot[7]>50) & (df_plot[7]<140)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_SS, [(rows[3],rows[2],rows[4],rows[6],rows[5],model)  for i,rows in df_plot.iterrows()])
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
        elif k==2:
           df_plot = df_plot[((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<165))]#['0','1','2','3','4',str(16+k*6+3),str(16+k*6+4),str(16+k*6+2)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_SSS, [(rows[3],rows[2],rows[4],rows[6],rows[5],model)  for i,rows in df_plot.iterrows()])
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
        elif k==3:
           df_plot = df_plot[((df_plot[7]>5) & (df_plot[7]<25)) | ((df_plot[7]>50) & (df_plot[7]<65)) ]#['0','1','2','3','4',str(16+k*6+3),str(16+k*6+4),str(16+k*6+2)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_S_ScS, [(rows[3],rows[2],rows[4],rows[6],rows[5],model,['ScS'])  for i,rows in df_plot.iterrows()])
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
        elif k==4:
           df_plot = df_plot[(df_plot[7]>100) & (df_plot[7]<150)]
           pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
           results=pool.starmap(points_Sdiff, [(rows[3],rows[2],rows[4],rows[6],rows[5],model)  for i,rows in df_plot.iterrows()])
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

        print(df_plot[["start_lon","start_lat","mid_lon","mid_lat"]])

        


        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
           
        
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+3],cmap=True,            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+3],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[16+k*6+3]),cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_amp_"+file_names[i]+"_"+ph+".png")
        plt.close()


        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
        #fig.coast(shorelines=True, frame=True)   
        
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+4],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.05c",fill=np.log(df_plot[16+k*6+4]),cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_env_"+file_names[i]+"_"+ph+".png")
        plt.close()



        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-10,10]) 
        fig.basemap(region="d",projection="N-150/12c")
        #fig.coast(shorelines=True, frame=True)   
        
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+2],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+2],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.05c",fill=df_plot[16+k*6+2],cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_time_"+file_names[i]+"_"+ph+".png")
        plt.close()

        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
        #fig.coast(shorelines=True, frame=True)   
        
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+2],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+2],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[16+k*6+5]),cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_amp_grl_"+file_names[i]+"_"+ph+".png")
        plt.close()

for i,dfs in enumerate(df):
    df_plot = dfs[(dfs[18]>max_cc[i]) & (dfs[23]>max_cc[i]) & (np.abs(dfs[19]<15)) & (np.abs(dfs[24]<15))]
    df_plot = df_plot[df_plot[8] == df_plot[8]]
    df_plot = df_plot[(df_plot[7]>50) & (df_plot[7]<140)]
    pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
    results=pool.starmap(points_SS, [(rows[3],rows[2],rows[4],rows[6],rows[5],model)  for i,rows in df_plot.iterrows()])
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
    fig=pygmt.Figure()
    pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
    fig.basemap(region="d",projection="N-150/12c")
    #fig.coast(shorelines=True, frame=True)   
    
    #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+3],cmap=True,            
    #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+3],cmap=True,transparency=40)
    fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[8]),cmap=True,transparency=40)
    fig.coast(shorelines=True, frame=True)
    fig.colorbar()
    fig.savefig(plot_dir+"/SS_S_"+file_names[i]+".png")
    plt.close()
    df_plot = dfs[(dfs[28]>max_cc[i]) & (dfs[23]>max_cc[i]) & (np.abs(dfs[29]<15)) & (np.abs(dfs[24]<15))]
    df_plot = df_plot[df_plot[10] == df_plot[10]]
    df_plot = df_plot[(((df_plot[7]>75) & (df_plot[7]<110)) | ((df_plot[7]>110) & (df_plot[7]<140)))]
    pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
    results=pool.starmap(points_SS, [(rows[3],rows[2],rows[4],rows[6],rows[5],model)  for i,rows in df_plot.iterrows()])
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
    fig=pygmt.Figure()
    pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
    fig.basemap(region="d",projection="N-150/12c")
    #fig.coast(shorelines=True, frame=True)   
    
    #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+3],cmap=True,            
    #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*6+3],cmap=True,transparency=40)
    fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[10]),cmap=True,transparency=40)
    fig.coast(shorelines=True, frame=True)
    fig.colorbar()
    fig.savefig(plot_dir+"/SSS_SS_"+file_names[i]+".png")
    plt.close()