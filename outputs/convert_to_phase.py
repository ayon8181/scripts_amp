import pandas as pd
import numpy as np


file_names=["real_data.txt","obsd_1D_crust.txt","obsd_3D_crust.txt","prem_3D_crust.txt","obsd_glad.txt"]#,"prem_3D.txt4","obsd_1D.txt4","obsd_3D.txt4"]
df=[]
for i,f in enumerate(file_names):
    temp    = pd.read_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/0_"+f,header=None,delimiter=" ",skipinitialspace=True)
    df.append(temp)

d_all=pd.merge(df[0],df[1],on=[0,1,2,3,4,5,6,7],suffixes=("_real","_1D"))
d_all=pd.merge(d_all,df[2],on=[0,1,2,3,4,5,6,7])
d_all=pd.merge(d_all,df[3],on=[0,1,2,3,4,5,6,7],suffixes=("_3D","_prem_3D"))
d_all=pd.merge(d_all,df[4],on=[0,1,2,3,4,5,6,7])
d_all=pd.merge(d_all,df[1],on=[0,1,2,3,4,5,6,7],suffixes=("_glad","_prem_atten"))
#d_all=pd.merge(d_all,df[6],on=[0,1,2,3,4,5,6,7])
#d_all=pd.merge(d_all,df[7],on=[0,1,2,3,4,5,6,7],suffixes=("_1D_atten","_3D_atten"))


phase_list=['S','SS','SSS','ScS','Sdiff']
names=["_real","_1D","_3D","_prem_3D","_glad"]#,"_prem_atten","_1D_atten","_3D_atten"]
for k,ph in enumerate(phase_list):
    col_list=[]
    for i in range(8):
        col_list.append(i)
    for n in names:
        for i in range(7):
            col_list.append(str(16+k*7+i)+n)
    #for i in range(6):
    #    col_list.append(16+k*6+i)
    print(d_all)
    df_out=d_all[col_list]
    df_out.to_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/"+ph+"_all.txt",index=False,na_rep=np.nan)

col_list=[]
for i in range(8):
    col_list.append(i)
for n in names:
    col_list.append("17"+n)
    col_list.append("18"+n)
    col_list.append("24"+n)
    col_list.append("25"+n)
    col_list.append("8"+n)
    col_list.append("9"+n)
"""
col_list.append(17)
col_list.append(18)
col_list.append(23)
col_list.append(24)
col_list.append(8)
col_list.append(9)   
"""
df_out=d_all[col_list]
df_out.to_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SS_S_all.txt",index=False,na_rep=np.nan)

col_list=[]
for i in range(8):
    col_list.append(i)
for n in names:
    col_list.append("31"+n)
    col_list.append("32"+n)
    col_list.append("24"+n)
    col_list.append("25"+n)
    col_list.append("10"+n)
    col_list.append("11"+n)
"""
col_list.append(29)
col_list.append(30)
col_list.append(23)
col_list.append(24)
col_list.append(10)
col_list.append(11)  
"""
df_out=d_all[col_list]
df_out.to_csv("/scratch1/09038/ayon8181/scripts_amp/outputs/SSS_SS_all.txt",index=False,na_rep=np.nan)