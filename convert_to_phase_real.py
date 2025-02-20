import pandas as pd
import numpy as np
from obspy.core.event import read_events
from geographiclib.geodesic import Geodesic
import csv
from obspy.imaging import beachball
import multiprocessing as mp
import os

ev_list=[]
with open("event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])
m_list=["872_","1482_","2405_","3326_","3941_","4534_","4922_","5405_","6026_","6868_","","9010_","11289_","12181_","13332_","15821_","15988_"]
n_list=[872,1482,2405,3326,3941,4534,4922,5405,6026,6868,8112,9010,11289,12181,13332,15821,15988]
file_names=[]
for j in m_list:
    file_names.append("20_40_obsd_3D_crust_"+j+"body_prem.txt")
def rad_pat_sh(rake,dip,az,st,tk_of):
    y = (np.cos(rake)*np.cos(dip)*np.cos(tk_of)*np.sin(az-st) +
       np.cos(rake)*np.sin(dip)*np.sin(tk_of)*np.cos(2*(az-st)) +
       np.sin(rake)*np.cos(2*dip)*np.cos(tk_of)*np.cos(az-st) -
       0.5*np.sin(rake)*np.sin(2*dip)*np.sin(tk_of)*np.sin(2*(az-st)))
    return y

def nodal_plane_filter(event):
    mt    = event.focal_mechanisms[0].moment_tensor.tensor
    #nod_p = beachball.mt2plane(beachball.MomentTensor(0,1,-1,0,0,0,26))
    nod_p = beachball.mt2plane(beachball.MomentTensor(mt.m_rr,mt.m_tt,mt.m_pp,mt.m_rt,mt.m_rp,mt.m_tp,0))
    strike= nod_p.strike
    dip   = nod_p.dip
    rake  = nod_p.rake
    l=0
    np_0=[]
    rps=[]
    for i in range(0,360):
        rp=np.abs(rad_pat_sh(np.deg2rad(rake),np.deg2rad(dip),np.deg2rad(i),np.deg2rad(strike),np.deg2rad(20.0)))
        rps.append(rp)
    rps=rps/max(rps)
    for i in range(0,360):
        if rps[i]<=0.1:
           np_0.append(i)       
    return np_0

def plane(nodal):
    if nodal<0:
       y=360+nodal
    elif nodal>=360:
         y=0+(nodal-360)
    else:
         y=nodal
    
    return y


def az_co(evla,evlo,stla,stlo,np_az):
    az=Geodesic.WGS84.Inverse(evla, evlo, stla, stlo)['azi1']
    #print(az,np_az)
    r=1
    for p in np_az:
            if plane(az) >= plane(p-10) and plane(az) <= plane(p+10):
                r=0
                #print(evla,evlo,stla,stlo)
    return r

def rad_pat_fltr(df):
    dfs=[]
    for events in ev_list:
        event=read_events("quakeml/"+events+".xml")[0]
        #print(event.origins)
        #print(event.focal_mechanisms)
        np=nodal_plane_filter(event)
        df_event=df[df[0]==events]
        fr=[]
        pool = mp.Pool(6)
        results=pool.starmap(az_co, [(rows[3],rows[2],rows[6],rows[5],np)  for i,rows in df_event.iterrows()])
        pool.close()
        pool.join()
        for r in results:
            if r is not None:
               fr.append(r)
        #print(fr)
        df_event = df_event.assign(r=pd.Series(fr).values)
        #print(df_event)
        df_event=df_event[df_event['r']==1]
        dfs.append(df_event)
    new_df=pd.DataFrame()
    
    new_df=pd.concat(list(dfs),ignore_index=True)
               
    return new_df


df=[]
for i,f in enumerate(file_names):
    temp    = pd.read_csv(f,header=None,delimiter=" ",skipinitialspace=True)
    df.append(temp)
df[0]=rad_pat_fltr(df[0])
df_glad=pd.read_csv("20_40_obsd_glad_body_prem.txt",header=None,delimiter=" ",skipinitialspace=True)
df_s40=pd.read_csv("20_40_obsd_3D_crust_body.txt",header=None,delimiter=" ",skipinitialspace=True)
df_real  =pd.read_csv("20_40_real_data_body_prem.txt",header=None,delimiter=" ",skipinitialspace=True)
# df_md2  =pd.read_csv("20_40_prem_modified_2_body.txt",header=None,delimiter=" ",skipinitialspace=True)
df_md3  =pd.read_csv("20_40_prem_modified_3_body.txt",header=None,delimiter=" ",skipinitialspace=True)
df_md3_2  =pd.read_csv("20_40_obsd_glad_body.txt",header=None,delimiter=" ",skipinitialspace=True) 
df_md = pd.read_csv("20_40_prem_modified_3_body_prem.txt",header=None,delimiter=" ",skipinitialspace=True)

#for i,f in enumerate(file_names):
#    temp    = pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/0_"+f,header=None,delimiter=" ",skipinitialspace=True)
#    df.append(temp)

d_all=pd.merge(df[0],df[1],how='left',on=[0,1],suffixes=("_"+str(n_list[0]),"_"+str(n_list[1])))

for i in range(2,len(m_list),2):    
    d_all=pd.merge(d_all,df[i],how='left',on=[0,1])
    if i+1 != len(m_list):
       d_all=pd.merge(d_all,df[i+1],how='left',on=[0,1],suffixes=("_"+str(n_list[i]),"_"+str(n_list[i+1])))
if i==len(m_list)-1:
    # d_all=pd.merge(d_all,df[-1],how='left',on=[0,1])
    d_all=pd.merge(d_all,df[-1],how='left',on=[0,1],suffixes=("_fow","_"+str(n_list[i])))
d_all=pd.merge(d_all,df_real,how='left',on=[0,1])
d_all=pd.merge(d_all,df_glad,how='left',on=[0,1],suffixes=("_real","_glad"))
d_all=pd.merge(d_all,df_md3,how='left',on=[0,1])
d_all=pd.merge(d_all,df_md3_2,how='left',on=[0,1],suffixes=("_mod_3","_m25"))
d_all=pd.merge(d_all,df_s40,how='left',on=[0,1])
d_all=pd.merge(d_all,df_md,how='left',on=[0,1],suffixes=("_s40","_mprem"))    
# d_all=pd.merge(d_all,df[4],how='left',on=[0,1])
# d_all=pd.merge(d_all,df[5],how='left',on=[0,1],suffixes=("_prem_16","_3D_atten"))
# d_all=pd.merge(d_all,df[6],how='left',on=[0,1])
# d_all=pd.merge(d_all,df[7],how='left',on=[0,1],suffixes=("_synt","_S80_3D"))
# d_all=pd.merge(d_all,df[8],how='left',on=[0,1])
# d_all=pd.merge(d_all,df[9],how='left',on=[0,1],suffixes=("_1D","_3D_11"))
# d_all=pd.merge(d_all,df[10],how='left',on=[0,1])
# d_all=pd.merge(d_all,df[11],how='left',on=[0,1],suffixes=("_3D_5","_3D_1"))
# d_all=pd.merge(d_all,df[12],how='left',on=[0,1])
# d_all=pd.merge(d_all,df[12],how='left',on=[0,1],suffixes=("_3D_18","_fow"))
  # No tag=prem_3d_atten
for i in range(2,8):
    d_all.rename(columns={str(i)+'_4922': i}, inplace=True)


phase_list=['S','SS','SSS','ScS']
names=[]
for j in n_list:
    names.append("_"+str(j))
for k,ph in enumerate(phase_list):
    col_list=[]
    for i in range(8):
        col_list.append(i)
    for n in names:
        for i in range(5):
            col_list.append(str(14+k*5+i)+n)
    for i in range(5):
        col_list.append(str(14+k*5+i)+"_glad")
        col_list.append(str(14+k*5+i)+"_real")
        #col_list.append(str(14+k*5+i)+"_mod_2")
        col_list.append(str(14+k*5+i)+"_mod_3")
        col_list.append(str(14+k*5+i)+"_m25")
        col_list.append(str(14+k*5+i)+"_s40")
        col_list.append(str(14+k*5+i)+"_mprem")
    #for i in range(6):
    #    col_list.append(16+k*6+i)
    print(d_all)
    df_out=d_all[col_list]
    df_out.to_csv("20_40_"+ph+"_all.txt",index=False,na_rep=np.nan)

col_list=[]
for i in range(8):
    col_list.append(i)
names.append("_glad")
names.append("_real")
#names.append("_mod_2")
names.append("_mod_3")
names.append("_m25")
names.append("_s40")
names.append("_mprem")
for n in names:
    
    col_list.append("8"+n)
    col_list.append("9"+n)
    col_list.append(str(14+0*5+2)+n)
    col_list.append(str(14+1*5+2)+n)
    col_list.append(str(14+1*5+4)+n)
    col_list.append(str(14+0*5+4)+n)
    col_list.append(str(14+1*5+1)+n)
    col_list.append(str(14+0*5+1)+n)
    col_list.append(str(14+1*5+0)+n)
    col_list.append(str(14+0*5+0)+n)
"""
col_list.append(17)
col_list.append(18)
col_list.append(23)
col_list.append(24)
col_list.append(8)
col_list.append(9)   
"""
df_out=d_all[col_list]
df_out.to_csv("20_40_SS_S_all.txt",index=False,na_rep=np.nan)

col_list=[]
for i in range(8):
    col_list.append(i)
for n in names:
    col_list.append("24"+n)
    col_list.append("25"+n)
    col_list.append("29"+n)
    col_list.append("30"+n)
    col_list.append("12"+n)
    col_list.append("13"+n)
"""
col_list.append(29)
col_list.append(30)
col_list.append(23)
col_list.append(24)
col_list.append(10)
col_list.append(11)  
"""
df_out=d_all[col_list]
df_out.to_csv("20_40_SSS_SS_all.txt",index=False,na_rep=np.nan)

# col_list=[]
# for i in range(8):
#     col_list.append(i)
# for n in names:
#     col_list.append("24"+n)
#     col_list.append("25"+n)
#     col_list.append("39"+n)
#     col_list.append("40"+n)
#     col_list.append("10"+n)
#     col_list.append("11"+n)
#     col_list.append(str(18+4*5+2)+n)
#     col_list.append(str(18+1*5+2)+n)
#     col_list.append(str(18+1*5+4)+n)
#     col_list.append(str(18+4*5+4)+n)
# """
# col_list.append(29)
# col_list.append(30)
# col_list.append(23)
# col_list.append(24)
# col_list.append(10)
# col_list.append(11)  
# """
# df_out=d_all[col_list]
# df_out.to_csv("20_40_SS_Sdiff_all.txt",index=False,na_rep=np.nan)


# col_list=[]
# for i in range(8):
#     col_list.append(i)
# for n in names:
#     col_list.append("29"+n)
#     col_list.append("30"+n)
#     col_list.append("39"+n)
#     col_list.append("40"+n)
#     col_list.append("16"+n)
#     col_list.append("17"+n)
# """
# col_list.append(29)
# col_list.append(30)
# col_list.append(23)
# col_list.append(24)
# col_list.append(10)
# col_list.append(11)  
# """
# df_out=d_all[col_list]
# df_out.to_csv("20_40_SSS_Sdiff_all.txt",index=False,na_rep=np.nan)