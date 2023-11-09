import pandas as pd
import numpy as np


file_names=["real_data","obsd_1D_crust","obsd_3D_crust","obsd_glad"]
df=[]
for i,f in enumerate(file_names):
    temp    = pd.read_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+f+".txt2",header=None,delimiter=" ",skipinitialspace=True)
    df.append(temp)

d_all=pd.merge(df[0],df[1],on=[0,1,2,3,4,5,6,7],suffixes=("_real","_1D"))
d_all=pd.merge(d_all,df[2],on=[0,1,2,3,4,5,6,7])
d_all=pd.merge(d_all,df[3],on=[0,1,2,3,4,5,6,7],suffixes=("_3D","_glad"))

phase_list=['S','SS','SSS','ScS','Sdiff']
names=["_real","_1D","_3D","_glad"]
for k,ph in enumerate(phase_list):
    col_list=[]
    for i in range(8):
        col_list.append(i)
    for n in names:
        for i in range(6):
            col_list.append(str(16+k*6+i)+n)
    print(d_all)
    df_out=d_all[col_list]
    df_out.to_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/"+ph+"_all.txt",index=False,na_rep=np.nan)

col_list=[]
for i in range(8):
    col_list.append(i)
for n in names:
    col_list.append("17"+n)
    col_list.append("23"+n)
    col_list.append("8"+n)
    col_list.append("9"+n)
df_out=d_all[col_list]
df_out.to_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SS_S_all.txt",index=False,na_rep=np.nan)

col_list=[]
for i in range(8):
    col_list.append(i)
for n in names:
    col_list.append("29"+n)
    col_list.append("23"+n)
    col_list.append("10"+n)
    col_list.append("11"+n)
df_out=d_all[col_list]
df_out.to_csv("/scratch1/09038/ayon8181/pypaw_workflow_test/outputs/SSS_SS_all.txt",index=False,na_rep=np.nan)