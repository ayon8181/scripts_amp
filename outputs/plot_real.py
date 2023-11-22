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
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 22

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)

ev_list=[]
with open("same.txt","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[2])
    
         

file_names      = ["_1D", "_3D","_glad","_real",""]
colors=['springgreen','red','blueviolet',"darkorange",'black']    
lss   =[2,2,2,2,2] 
style =['None','solid',"dashed","dotted","dashdot"]
labels=['S40RTS_1D_Crust', 'S40RTS_3D_crust','GLAD_M25',"real_data","prem_3D_crust"] 
fccs  =['springgreen','None','None','None',"None"] 
df=[]
plot_dir="/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/plots_map_real"
model=TauPyModel(model="prem")
#for i,f in enumerate(file_names):
#    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt",header=None,delimiter=" ",skipinitialspace=True)
 #   df.append(temp)
def points_S_ScS(evla, evlo, evdp, stla, stlo, model, phase_list):
    #print(evla,evlo,stla,stlo)
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
    df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+ph+"_all.txt",skipinitialspace=True,delimiter=",")
    df_plot=df[df[str(16+k*6+5)+"_real"].notna()]
    df_plot=df_plot[(np.abs(df_plot[str(16+k*6+2)+"_real"])<15) & (np.abs(df_plot[str(16+k*6+2)+"_1D"])<15) &(np.abs(df_plot[str(16+k*6+2)+"_3D"])<15) & (np.abs(df_plot[str(16+k*6+2)+"_glad"])<15)]
    df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    #df_plot=df_plot[df_plot["4"]>100]
    
       # df_plot = dfs[(dfs[16+k*5+1]>max_cc[i]) & (np.abs((dfs[16+k*5+2])<15))]
        #df_plot = df_plot[df_plot[16+k*5+3] == df_plot[16+k*5+3]]
    if k==0:
        df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points_S_ScS, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],model,['S'])  for i,rows in df_plot.iterrows()])
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
        df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points_SS, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],model)  for i,rows in df_plot.iterrows()])
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
        df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<165))]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points_SSS, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],model)  for i,rows in df_plot.iterrows()])
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
        df_plot = df_plot[((df_plot["7"]>5) & (df_plot["7"]<25)) | ((df_plot["7"]>50) & (df_plot["7"]<65)) ]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points_S_ScS, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],model,['ScS'])  for i,rows in df_plot.iterrows()])
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
        df_plot = df_plot[(df_plot["7"]>100) & (df_plot["7"]<150)]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points_Sdiff, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],model)  for i,rows in df_plot.iterrows()])
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

        
    for i,df in enumerate(file_names):

        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
        
        
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(16+k*6+5)+df]),cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_amp_"+file_names[i]+"_"+ph+".png")
        plt.close()


        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
        #fig.coast(shorelines=True, frame=True)   
        
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(16+k*6+4)+df]),cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_env_"+file_names[i]+"_"+ph+".png")
        plt.close()



        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c")
        #fig.coast(shorelines=True, frame=True)   
        
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+2],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+2],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(16+k*6+3)+df]),cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        fig.colorbar()
        fig.savefig(plot_dir+"/global_amp_3d_1d_"+file_names[i]+"_"+ph+".png")
        plt.close()
    plt.figure(1,figsize=(10,10))
    for i,df in enumerate(file_names):
        plt.hist(np.log(df_plot[str(16+k*6+5)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        plt.axvline(np.mean(np.log(df_plot[str(16+k*6+5)+df])),color=colors[i],linestyle="dashed",linewidth=2)
    plt.legend()
    plt.xlabel("Log of Amp Ratio Measurements for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_amp_r_"+ph+".png")
    plt.close()
    for i,df in enumerate(file_names):
        plt.hist(np.log(df_plot[str(16+k*6+3)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        plt.axvline(np.mean(np.log(df_plot[str(16+k*6+3)+df])),color=colors[i],linestyle="dashed",linewidth=2)
    plt.legend()
    plt.xlabel("Log of 3d/1d amplitude ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_3d_1d_"+ph+".png")
    plt.close()
    for i,df in enumerate(file_names):
        plt.hist(np.log(df_plot[str(16+k*6+4)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
        plt.axvline(np.mean(np.log(df_plot[str(16+k*6+4)+df])),color=colors[i],linestyle="dashed",linewidth=2)
    plt.legend()
    plt.xlabel("Log of 3d/1d envp ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_envp_3d_1d_"+ph+".png")
    plt.close()


df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_plot=df[df[str(8)+"_real"].notna()]
df_plot = df_plot[(df_plot["7"]>50) & (df_plot["7"]<140)]
df_plot=df_plot[(np.abs(df_plot[str(18)+"_real"])<15) & (np.abs(df_plot[str(18)+"_1D"])<15) &(np.abs(df_plot[str(18)+"_3D"])<15) & (np.abs(df_plot[str(18)+"_glad"])<15)]
df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<15) & (np.abs(df_plot[str(24)+"_1D"])<15) &(np.abs(df_plot[str(24)+"_3D"])<15) & (np.abs(df_plot[str(24)+"_glad"])<15)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(points_SS, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],model)  for i,rows in df_plot.iterrows()])
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
for i,f in enumerate(file_names):
    fig=pygmt.Figure()
    pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
    fig.basemap(region="d",projection="N-150/12c")
    #fig.coast(shorelines=True, frame=True)   

    #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
    #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
    fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(8)+f]),cmap=True,transparency=40)
    fig.coast(shorelines=True, frame=True)
    fig.colorbar()
    fig.savefig(plot_dir+"/SS_S_"+file_names[i]+".png")
    plt.close()
plt.figure(1,figsize=(10,10))
for i,df in enumerate(file_names):
    plt.hist(np.log(df_plot[str(8)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
    plt.axvline(np.mean(np.log(df_plot[str(8)+df])),color=colors[i],linestyle="dashed",linewidth=2)
plt.legend()
plt.xlabel("Log of SS/S Measurements ")
plt.ylabel("Counts")
plt.savefig(plot_dir+"/Histogram_SS_S.png")
plt.close()

df=pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SSS_SS_all.txt",skipinitialspace=True,delimiter=",")
df_plot=df[df[str(10)+"_real"].notna()]
df_plot = df_plot[(((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<140)))]
df_plot=df_plot[(np.abs(df_plot[str(30)+"_real"])<15) & (np.abs(df_plot[str(30)+"_1D"])<15) &(np.abs(df_plot[str(30)+"_3D"])<15) & (np.abs(df_plot[str(30)+"_glad"])<15)]
df_plot=df_plot[(np.abs(df_plot[str(24)+"_real"])<15) & (np.abs(df_plot[str(24)+"_1D"])<15) &(np.abs(df_plot[str(24)+"_3D"])<15) & (np.abs(df_plot[str(24)+"_glad"])<15)]
pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
results=pool.starmap(points_SS, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],model)  for i,rows in df_plot.iterrows()])
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
for i,f in enumerate(file_names):
    fig=pygmt.Figure()
    pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
    fig.basemap(region="d",projection="N-150/12c")
    #fig.coast(shorelines=True, frame=True)   
    
    #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
    #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.001p,+z",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
    fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(10)+f]),cmap=True,transparency=40)
    fig.coast(shorelines=True, frame=True)
    fig.colorbar()
    fig.savefig(plot_dir+"/SSS_SS_"+file_names[i]+".png")
    plt.close()

plt.figure(1,figsize=(10,10))
for i,df in enumerate(file_names):
    plt.hist(np.log(df_plot[str(10)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = 0.5)
    plt.axvline(np.mean(np.log(df_plot[str(10)+df])),color=colors[i],linestyle="dashed",linewidth=2)
plt.legend()
plt.xlabel("Log of SSS/SS Measurements")
plt.ylabel("Counts")
plt.savefig(plot_dir+"/Histogram_SSS_SS.png")
plt.close()