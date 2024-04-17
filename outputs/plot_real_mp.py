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
#import seaborn as sns
import scipy.stats as stats
import geopy.distance
from geographiclib.geodesic import Geodesic
from matplotlib.lines import Line2D
#from mpi4py import MPI

font = { "family":"serif",
          "color": "darkred",
          "weight":"normal",
          }

SMALL_SIZE = 16
MEDIUM_SIZE = 24
BIGGER_SIZE = 50

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

ev_list=[]
with open("./same.txt","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])


    
         

file_names      = ["_real","_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['darkred','cyan','yellowgreen','chocolate','black']    
lss   =[2,2.5,2,1.5,1] 
style =['None','solid',"solid","solid","solid"]
labels=["real_data",'S40RTS_1D_Crust','S40RTS_3D_crust',"GLAD_M25","PREM_Crust2.0"] #"PREM_QRFSI12",
fccs  =['hotpink','None','None','None',"None"] 
markers=[".","+","o","*","x"]
alp   =[0.5,1,1,1,1]
df=[]
plot_dir="/scratch1/09038/ayon8181/scripts_amp/outputs/plots_map/same"
model=TauPyModel(model="prem")
#for i,f in enumerate(file_names):
#    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt",header=None,delimiter=" ",skipinitialspace=True)
 #   df.append(temp)
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

handles=[]
for i,df in enumerate(file_names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker=markers[i], label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends")
plt.close()


phase_list=['S','SS','SSS','ScS','Sdiff']
for k,ph in enumerate(phase_list):
    df=pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all.txt",skipinitialspace=True,delimiter=",")
    df_plot=df[df[str(16+k*7+5)+"_real"].notna()]
    df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    #df_plot=df_plot[(np.abs(df_plot[str(16+k*7+2)+"_real"])<15) & (np.abs(df_plot[str(16+k*7+2)+"_1D"])<15) &(np.abs(df_plot[str(16+k*7+2)+"_3D"])<15) & (np.abs(df_plot[str(16+k*7+2)+"_prem_3D"])<15)]
    

    
    if k==0:    
        df_plot = df_plot[(df_plot["7"]>30) & (df_plot["7"]<70)]
        print(df_plot)#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],2/5.0,3/5.0)  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        """
        comm =MPI.COMM_WORLD
        size =comm.Get_size()
        rank =comm.Get_rank()

        inputs=[(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],2/5.0,3/5.0)  for i,rows in df_plot.iterrows()]

        inputs_chunk=comm.scatter(inputs,root=0)

        results_chunk=[points(*args) for args in inputs_chunk]

        results = comm.gather(results_chunk, root=0)

        if rank==0:
            results = [result for chunk in results for result in chunk]
        
        MPI.Finalize()
        """


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
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        """
        comm =MPI.COMM_WORLD
        size =comm.Get_size()
        rank =comm.Get_rank()

        inputs=[(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()]

        inputs_chunk=comm.scatter(inputs,root=0)

        results_chunk=[points(*args) for args in inputs_chunk]

        results = comm.gather(results_chunk, root=0)

        if rank==0:
            results = [result for chunk in results for result in chunk]
        MPI.Finalize()

        """
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
        df_plot = df_plot[((df_plot["7"]>75) & (df_plot["7"]<110)) | ((df_plot["7"]>110) & (df_plot["7"]<150))]#['0','1','2','3','4',str(16+k*5+3),str(16+k*5+4),str(16+k*5+2)]
        pool = mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK']))
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        """
        comm =MPI.COMM_WORLD
        size =comm.Get_size()
        rank =comm.Get_rank()

        inputs=[(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()]

        inputs_chunk=comm.scatter(inputs,root=0)

        results_chunk=[points(*args) for args in inputs_chunk]

        results = comm.gather(results_chunk, root=0)

        if rank==0:
            results = [result for chunk in results for result in chunk]
        MPI.Finalize()
        """

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
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        """
        comm =MPI.COMM_WORLD
        size =comm.Get_size()
        rank =comm.Get_rank()

        inputs=[(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()]

        inputs_chunk=comm.scatter(inputs,root=0)

        results_chunk=[points(*args) for args in inputs_chunk]

        results = comm.gather(results_chunk, root=0)

        if rank==0:
            results = [result for chunk in results for result in chunk]
        MPI.Finalize()
        """
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
        results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()])
        pool.close()
        pool.join()
        """
        comm =MPI.COMM_WORLD
        size =comm.Get_size()
        rank =comm.Get_rank()

        inputs=[(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],45/100.0,55/100.0)  for i,rows in df_plot.iterrows()]

        inputs_chunk=comm.scatter(inputs,root=0)

        results_chunk=[points(*args) for args in inputs_chunk]

        results = comm.gather(results_chunk, root=0)

        if rank==0:
            results = [result for chunk in results for result in chunk]
        MPI.Finalize()
        """
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

    

       
    for i,df in enumerate(file_names):

        print(df_plot)

        ### Plotting Global Maps for Amplitude Ratio
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        
        
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(16+k*7+5)+df]),pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        if i==0:
           fig.colorbar()
        fig.savefig(plot_dir+"/global_amp_"+file_names[i]+"_"+ph+".png",dpi=600)
        plt.close()

        ### Plotting Global Maps for Envelope Ratio
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        #fig.coast(shorelines=True, frame=True)   
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=df_plot[str(16+k*7+4)+df],pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        if i==0:
           fig.colorbar()
        fig.savefig(plot_dir+"/global_env_"+file_names[i]+"_"+ph+".png")
        plt.close()
        
        ### Plotting Global Maps for CC Traveltime
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        #fig.coast(shorelines=True, frame=True)   
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=np.log(df_plot[str(16+k*7+3)+df]),pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        if i==0:
           fig.colorbar()
        fig.savefig(plot_dir+"/global_3d1d_"+file_names[i]+"_"+ph+".png")
        plt.close()
        
        ### Plotting Global Maps for 3D/1D Amplitude
        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-10,10]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        #fig.coast(shorelines=True, frame=True)   
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+2],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+2],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=(df_plot[str(16+k*7+2)+df]),pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        if i==0:
           fig.colorbar()
        fig.savefig(plot_dir+"/global_cc_t_"+file_names[i]+"_"+ph+".png")
        plt.close()

        fig=pygmt.Figure()
        pygmt.makecpt(cmap='polar',reverse=True,series=[-1.5,1.5]) 
        fig.basemap(region="d",projection="N-150/12c",frame=True)
        #fig.coast(shorelines=True, frame=True)   
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)#zvalue=df_plot[16+k*5+3],cmap=True,            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000005p",transparency=40)
        #fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
        #fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.000001p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.1c",fill=df_plot[str(16+k*7+6)+df],pen="0.000001p",cmap=True,transparency=40)
        fig.coast(shorelines=True, frame=True)
        if i==0:
           fig.colorbar()
        fig.savefig(plot_dir+"/global_amp_misfit_"+file_names[i]+"_"+ph+".png")
        plt.close()
        

    #sns.set_theme(style="darkgrid")
    #fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
    #df_sns=pd.DataFrame()
        

    ## Plot scatter plots for epicentral distance and amplitude ratios



    ## Plot Histograms for Amplitude Ratios
    plt.figure(1,figsize=(7.08,7.08))
    for i,df in enumerate(file_names):
        
       
        plt.hist(np.log(df_plot[str(16+k*7+5)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i], linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])

    plt.xlim(-1.5,1.5)
    plt.legend(fontsize=12)
    plt.xlabel("$X_{2}$ for "+ph+" phase")
    plt.ylabel("Counts")
    plt.tight_layout()
    plt.savefig(plot_dir+"/Histogram_amp_r_"+ph+".png",dpi=600)
    plt.close()
    

    ### Plot Histograms for 3D/1D Amplitude Ratios
    plt.figure(1,figsize=(7.08,7.08))
    for i,df in enumerate(file_names):
        

        plt.hist(np.log(df_plot[str(16+k*7+3)+df]), bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=lss[i], linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])

    plt.xlim(-1.5,1.5)
    plt.legend(fontsize=12)
    plt.xlabel("Energy Ratio for "+ph+" phase")
    plt.ylabel("Counts")
    plt.tight_layout()
    plt.savefig(plot_dir+"/Histogram_amp_3d1d_"+ph+".png",dpi=600)
    plt.close()
    
    # Plot Histograms for Cross Correlation Traveltime
    plt.figure(1,figsize=(7.08,7.08))
    for i,df in enumerate(file_names):
        plt.hist(df_plot[str(16+k*7+2)+df], bins=23, range=(-10,10), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])
        #bin_centers = 0.5 * (bins[:-1] + bins[1:])
        #plt.plot(bin_centers, values, color=colors[i],linewidth=1)
        #plt.axvline(np.mean(np.log(df_plot[str(16+k*7+3)+df])),color=colors[i],linewidth=2)
    plt.xlim(-10,10)
    plt.legend(fontsize=12)
    plt.xlabel("Cross-correlation Travel Time Shift for "+ph+" phase")
    plt.ylabel("Counts")
    plt.tight_layout()
    plt.savefig(plot_dir+"/Histogram_cc_t_"+ph+".png",dpi=600)
    plt.close()

    # Plot Histograms for Envelope Ratios
    plt.figure(1,figsize=(7.08,7.08))
    for i,df in enumerate(file_names):
        plt.hist(df_plot[str(16+k*7+4)+df], bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])
        #bin_centers = 0.5 * (bins[:-1] + bins[1:])
        #plt.plot(bin_centers, values, color=colors[i],linewidth=1)
        #plt.axvline(np.mean(np.log(df_plot[str(16+k*7+4)+df])),color=colors[i],linewidth=2)
    plt.xlim(-1.5,1.5)
    plt.legend(fontsize=12)
    plt.xlabel("$X_{env}$ for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_envp_3d_1d_"+ph+".png",dpi=600)
    plt.close()

    plt.figure(1,figsize=(7.08,7.08))
    for i,df in enumerate(file_names):
        plt.hist(df_plot[str(16+k*7+6)+df], bins=23, range=(-1.5,1.5), edgecolor=colors[i], linewidth=3, linestyle=style[i],label=labels[i],facecolor=fccs[i], alpha = alp[i])
        #bin_centers = 0.5 * (bins[:-1] + bins[1:])
        #plt.plot(bin_centers, values, color=colors[i],linewidth=1)
        #plt.axvline(np.mean(np.log(df_plot[str(16+k*7+4)+df])),color=colors[i],linewidth=2)
    plt.xlim(-1.5,1.5)
    plt.legend(fontsize=12)
    plt.xlabel("$X_{1}$ for "+ph+" phase")
    plt.ylabel("Counts")
    plt.savefig(plot_dir+"/Histogram_amp_misfit_"+ph+".png",dpi=600)
    plt.close()


    ## Plot Line Histograms
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
       
        y, edges = np.histogram(np.log(df_plot[str(16+k*7+3)+df]), bins=23, range=(-1.5,1.5))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i])
    plt.text(1.3,0.8*max(y/np.sum(y)),ph)
    plt.xlim(-1.5,1.5)
    #plt.legend(fontsize=12)
    plt.xlabel("$X_{2}$ for "+ph+" phase")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/line_amp_r_"+ph+".png",dpi=600)
    plt.close()
    

    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        
        y, edges = np.histogram(np.log(df_plot[str(16+k*7+3)+df]), bins=23, range=(-1.5,1.5))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i])

    plt.text(1.3,0.8*max(y/np.sum(y)),ph)
    plt.xlim(-1.5,1.5)
    #plt.legend(fontsize=12)
    plt.xlabel("Energy Ratio for "+ph+" phase")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/line_amp_3d1d_"+ph+".png",dpi=600)
    plt.close()
    
    # Plot Histograms for Cross Correlation Traveltime
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        y, edges = np.histogram(df_plot[str(16+k*7+2)+df], bins=23, range=(-10,10))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i])
    plt.text(8.5,0.8*max(y/np.sum(y)),ph)
    plt.xlim(-10,10)
    #plt.legend(fontsize=12)
    plt.xlabel("Cross-correlation Travel Time Shift for "+ph+" phase")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/line_cc_t_"+ph+".png",dpi=600)
    plt.close()

    # Plot Histograms for Envelope Ratios
    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        y, edges = np.histogram(df_plot[str(16+k*7+4)+df], bins=23, range=(-1.5,1.5))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i])
    plt.text(1.3,0.8*max(y/np.sum(y)),ph)
    plt.xlim(-1.5,1.5)
    #plt.legend(fontsize=12)
    plt.xlabel("$X_{env}$ for "+ph+" phase")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/line_envp_3d_1d_"+ph+".png",dpi=600)
    plt.close()

    plt.figure(1,figsize=(7.08,3.54))
    for i,df in enumerate(file_names):
        y, edges = np.histogram(df_plot[str(16+k*7+6)+df], bins=23, range=(-1.5,1.5))
        centers = 0.5 * (edges[1:] + edges[:-1])
        plt.plot(centers,y/np.sum(y),marker=markers[i],color=colors[i],alpha=alp[i],lw=lss[i])
    plt.text(1.3,0.8*max(y/np.sum(y)),ph)
    plt.xlim(-1.5,1.5)
    #plt.legend(fontsize=12)
    plt.xlabel("$X_{1}$ for "+ph+" phase")
    plt.ylabel("Normalized Counts")
    plt.tight_layout()
    plt.grid()
    plt.savefig(plot_dir+"/line_amp_misfit_"+ph+".png",dpi=600)
    plt.close()
    

    for i,df in enumerate(file_names[1:]):
       
       fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
       
       ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.scatter(df_plot[str(16+k*7+4)+"_real"],df_plot[str(16+k*7+4)+df],marker=".",alpha=0.6)
       ax.set_xlabel("$X_{env}$ for Real Data")
       ax.set_ylabel("$X_{env}$ for "+labels[i+1])
       
       ax.set_xlim([-1.5, 1.5])
       ax.set_ylim([-1.5, 1.5])
       plt.tight_layout()
       plt.savefig(plot_dir+"/Scatter_envp_3d_1d_"+ph+"_"+labels[i+1]+".png",dpi=600)
       plt.close()
       
       
       #slope, intercept, r, p, se = stats.linregress(x=np.log(df_plot[str(16+k*7+5)+"_real"]), y=np.log(df_plot[str(16+k*7+5)+df]))
       fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
       #sns.regplot(x=np.log(df_plot[str(16+k*7+5)+"_real"]), y=np.log(df_plot[str(16+k*7+5)+df]),  color="m",ax=ax)
       ax.set_xlim(-1.5,1.5)
       ax.set_ylim(-1.5,1.5)
       ax.set_xlabel("$X_{2}$ for Real Data "+ph,fontdict=font)
       ax.set_ylabel("$X_{2}$ for "+labels[i+1]+" "+ph,fontdict=font)
       ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.scatter(np.log(df_plot[str(16+k*7+5)+"_real"]),np.log(df_plot[str(16+k*7+5)+df]),marker=".",alpha=0.6)
       plt.tight_layout()
       plt.savefig(plot_dir+"/Scatter_amp_r_"+ph+"_"+labels[i+1]+".png",dpi=600)
       plt.close()

       
       fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
       #ax1=sns.kdeplot(x=np.log(df_plot[str(16+k*7+3)+"_real"]), y=np.log(df_plot[str(16+k*7+3)+df]), cmap='viridis', clip=(-1.5, 1.5), shade=True)
       ax.plot(np.linspace(-10,10,100), np.linspace(-10,10,100), linestyle='--', color='black')
       ax.plot(np.linspace(-10,10,100), np.linspace(-10,10,100), linestyle='--', color='black')
       ax.scatter(df_plot[str(16+k*7+2)+"_real"],df_plot[str(16+k*7+2)+df],marker=".",alpha=0.6)
       ax.set_xlabel("Cross Correlation Traveltime for Real Data")
       ax.set_ylabel("Cross Correlation Traveltime for "+labels[i+1])
       ax.set_xlim([-10, 10])
       ax.set_ylim([-10, 10])
       plt.tight_layout()
       plt.savefig(plot_dir+"/Scatter_cc_t_"+ph+"_"+labels[i+1]+".png",dpi=600)
      
       plt.close()

       fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
       
       ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.scatter(df_plot[str(16+k*7+6)+"_real"],df_plot[str(16+k*7+3)+df],marker=".",alpha=0.6)
       ax.set_xlabel("$X_{1}$ for Real Data")
       ax.set_ylabel("$X_{1}$ for "+labels[i+1])
       
       ax.set_xlim([-1.5, 1.5])
       ax.set_ylim([-1.5, 1.5])
       plt.tight_layout()
       plt.savefig(plot_dir+"/Scatter_amp_misfit_"+ph+"_"+labels[i+1]+".png",dpi=600)
       plt.close()


       fig,ax=plt.subplots(1,1,figsize=(7.08,7.08))
       #sns.regplot(x=np.log(df_plot[str(16+k*7+5)+"_real"]), y=np.log(df_plot[str(16+k*7+5)+df]),  color="m",ax=ax)
       ax.set_xlim(-1.5,1.5)
       ax.set_ylim(-1.5,1.5)
       ax.set_xlabel("Energy Ratio for Real Data "+ph,fontdict=font)
       ax.set_ylabel("Energy Ratio for "+labels[i+1]+" "+ph,fontdict=font)
       ax.plot(np.linspace(-1.5,1.5,100), np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.plot(np.linspace(-1.5,1.5,100), -np.linspace(-1.5,1.5,100), linestyle='--', color='black')
       ax.scatter(np.log(df_plot[str(16+k*7+3)+"_real"]),np.log(df_plot[str(16+k*7+3)+df]),marker=".",alpha=0.6)
       plt.tight_layout()
       plt.savefig(plot_dir+"/Scatter_amp_3d_1d_"+ph+"_"+labels[i+1]+".png",dpi=600)
       plt.close()
  
