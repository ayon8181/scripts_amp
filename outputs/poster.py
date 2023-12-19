from pyasdf import ASDFDataSet
import obspy
import matplotlib.pyplot as plt
from scipy import signal

labels=['S40RTS_1D_Crust', 'S40RTS_3D_crust',"GLAD_M25","real_data"]
names=["proc_obsd_1D_crust","proc_obsd_3D_crust","proc_obsd_glad","proc_real_data"]#,"proc_prem_3D_crust","proc_obsd_glad","proc_real_data"]
tags =["proc_obsd_1","proc_obsd_3","proc_obsd_25","proc_real_data"]
colors=["brown","cyan","yellowgreen","blueviolet"]

S=[[0.25,1.25,0.25,1.75],[-3.75,-2.5,-4.75,-1.0],[-4.25,-0.75,-7.0,-8.25]]
arr=[1374,1704,1917]
ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/200904060132A.T017-040s.proc_synt.h5",mode="r")
data=ds.waveforms['TA.332A'].proc_synt.select(component="T")[0]
sampling_rate = data.stats.sampling_rate

st_ind=int(1300*sampling_rate)
en_ind=int(2000*sampling_rate)
p_data=data[st_ind:en_ind]
time=data.times()[st_ind:en_ind]

plt.figure(1,figsize=(15,5))
plt.plot(time,p_data,color="black",lw=1,label="1D PREM")
for i in range(4):
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/200904060132A.T017-040s."+names[i]+".h5",mode="r")
    data=ds.waveforms['TA.332A'][tags[i]].select(component="T")[0]
    p_data=data[st_ind:en_ind]
    taper=signal.windows.tukey(len(p_data),alpha=0.05)
    plt.plot(time,p_data*taper,color=colors[i],lw=1,label=labels[i])
plt.legend(loc="upper right")
plt.xlabel("Time (s)")
plt.savefig("poster_seis.svg",dpi=600)
plt.close()

fig,ax=plt.subplots(1,3,figsize=(30,5))
for j in range(3):
    ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/200904060132A.T017-040s.proc_synt.h5",mode="r")
    data=ds.waveforms['TA.332A'].proc_synt.select(component="T")[0]
    sampling_rate = data.stats.sampling_rate

    st_ind=int((arr[j]-50)*sampling_rate)
    en_ind=int((arr[j]+50)*sampling_rate)
    p_data=data[st_ind:en_ind]
    taper=signal.windows.tukey(len(p_data),alpha=0.05)
    time=data.times()[st_ind:en_ind]

   
    ax[j].plot(time,p_data*taper,color="black",lw=1,label="1D PREM")
    for i in range(4):
        ds=ASDFDataSet("/scratch1/09038/ayon8181/pypaw_workflow_test/seis/proc/200904060132A.T017-040s."+names[i]+".h5",mode="r")
        data=ds.waveforms['TA.332A'][tags[i]].select(component="T")[0]
        p_data=data[st_ind+int(S[j][i]*sampling_rate):en_ind+int(S[j][i]*sampling_rate)]
        taper=signal.windows.tukey(len(p_data),alpha=0.05)
        ax[j].plot(time,p_data*taper,color=colors[i],lw=1,label=labels[i])
    ax[j].legend(loc="upper right")
    ax[j].set_xlabel("Time (s)")
plt.savefig("poster_seis_2.svg",dpi=600)
plt.close()


