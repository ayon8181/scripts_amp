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
import matplotlib.gridspec as gridspec
#import seaborn as sns
import scipy.stats as stats
import geopy.distance
from geographiclib.geodesic import Geodesic
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
#from mpi4py import MPI

font = { "family":"serif",
          "color": "darkred",
          "weight":"normal",
          }

SMALL_SIZE = 28
MEDIUM_SIZE = 28
BIGGER_SIZE = 50

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)

font = {'family' : 'serif',
         'serif':  'cmr10'
         }

plt.rc('font', **font,size=SMALL_SIZE)  
plt.rc('axes.formatter', use_mathtext=True)


plt.rcParams['mathtext.default'] = 'rm'
plt.rc("mathtext", fontset="cm")
ev_list=[]
with open("event_list","r") as txt:
     reader=csv.reader(txt,skipinitialspace=True,delimiter="/")
     for row in reader:
         ev_list.append(row[0])




us=["NN","CI","TX","UW","LD","NM","N4","TA","WU","NC","AK","ET","UU"]
gsn=["II","IU","CU","IC","GT","US","CN"]         

file_names      = ["_real","_1D","_3D","_3D_2","_3D_15","_prem_3D","_glad","_synt",]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_prem_atten","_1D_atten","_3D_atten"]#"_1D","_3D","_glad","_prem_3D"]#"_1D","_3D","_glad","_prem_3D"]#
colors=['olive','crimson','dodgerblue','darkgreen','black','peru','indianred','grey','black']    
lss   =[1,1,1,1,1,1,1,1]
style =['None','solid',"solid","solid","solid","solid",'solid']
labels=["real_data",'SPREMc','S8000',"S11000","S5000","PREMc","GLAD_M25"] #"PREM_QRFSI12",
fccs  =['None','None','None','None',"None","None"] 
markers=[".","+","o","*","x","^","v","8",">","<","6"]
alp   =[0.5,1,1,1,1,1,1]
df=[]
plot_dir="new_figures_test"
model=TauPyModel(model="prem")

file_names=[]
labels=[]
n_list=[872,1482,2405,3326,3941,4534,4922,5405,6026,8112,11289,12181,13332,15821,15988]
for n in n_list:
    labels.append(str(n))
    file_names.append("_"+str(n))
file_names.append("_glad")
file_names.append("_real")
#file_names.append("_mod")
file_names.append("_mod_3")

def exp_format(y, _):
    return f"{np.exp(y):.2f}"


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
"""
handles=[]
for i,df in enumerate(file_names):
    handles.append(Line2D([0], [0], color=colors[i], lw=lss[i], marker=".", label=labels[i]))
                                        
fig, ax = plt.subplots()
ax.set_axis_off()
ax.legend(handles=handles, loc='center')
plt.savefig(plot_dir+"/legends")
plt.close()
"""

phase_list=['SS_S','SS_Sdiff']
ks=[0,1,2,4]
o_l=[8,10,12,34,16]
file_names2=["_real","_8112","_mod_3"]
labels=["Observed","S40RTS","mPREM"]


df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
    #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
df_plot=df[(df["0"].isin(ev_list))]

df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<140)]
df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
~(df_plot["1"].str.split(".").str[0].isin(us)) &
~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
(df_plot["1"].str.split(".").str[0].isin(gsn))]
df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
for i,df in enumerate(file_names):
    for w in [0,1]:
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
        df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]
    
fig,ax=plt.subplots(1,2,figsize=(15,7))
fig.text(0.02, 0.92, 'a)', fontsize=28, fontweight='bold')  # Position of "a)"
#fig.text(0.02, 0.5, 'c)', fontsize=28, fontweight='bold') 
fig.text(0.5, 0.92, 'b)', fontsize=28, fontweight='bold')
#fig.text(0.5, 0.5, 'd)', fontsize=28, fontweight='bold')

for i,f in enumerate(file_names2):
    ax[0].scatter(df_plot["7"],df_plot[str(14+1*5+2)+f]-df_plot[str(14+0*5+2)+f],alpha=0.5,marker=".",color="blue")
    ax[1].scatter(df_plot["7"],df_plot[str(14+1*5+4)+f]-df_plot[str(14+0*5+4)+f],alpha=0.5,marker=".",color="blue")
    fig.supxlabel(r"Epicentral Distaces$\degree$")
    ax[0].set_ylabel(r"$\delta \rm T^{ss-s}$")
    ax[1].set_ylabel(r"$\delta \rm lnA^{ss-s}$")
    ax[0].set_xlim(55,140)
    ax[1].set_xlim(55,140)
    ax[0].set_ylim(-12,12)
    ax[1].set_ylim(-2,2)
    ax[0].set_xticks(np.arange(55,141,10))
    ax[1].set_xticks(np.arange(55,141,10))
    ax[0].set_yticks(np.arange(-12,13,4))
    ax[1].set_yticks(np.arange(-2,2.1,0.8))
    ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(2))
    ax[1].yaxis.set_minor_locator(ticker.MultipleLocator(0.2))

fig.tight_layout()
fig.savefig(plot_dir+"/relative/scatters/SS_S_distance.pdf",dpi=600)
plt.close(fig)
fig = plt.figure(figsize=(14, 15))  
gs = gridspec.GridSpec(5, 6, height_ratios=[1, 0.3, 1, 0.05, 1], hspace=0.1, wspace=0.4)   
ax = [0, 0, 0, 0, 0,0,0]

# ax[6] = fig.add_subplot(gs[0, :])  # Top row spanning both columns
# ax[0] = fig.add_subplot(gs[2, 0])  # First subplot in second row (after space)
# ax[1] = fig.add_subplot(gs[2, 1])  # Second subplot in second row
# ax[2] = fig.add_subplot(gs[2, 2])  # First subplot in third row
# ax[3] = fig.add_subplot(gs[4, 0])  # First subplot in third row
# ax[4] = fig.add_subplot(gs[4, 1]) 
# ax[5] = fig.add_subplot(gs[4, 2])
ax[6] = fig.add_subplot(gs[0, :])  # Top row spanning both columns
ax[0] = fig.add_subplot(gs[2, 0:2])  # First subplot in second row (after space)
ax[1] = fig.add_subplot(gs[2, 2:4])  # Second subplot in second row
ax[2] = fig.add_subplot(gs[2, 4:6])  # First subplot in third row
ax[3] = fig.add_subplot(gs[4, 0:2])  # First subplot in third row
ax[4] = fig.add_subplot(gs[4, 2:4]) 
ax[5] = fig.add_subplot(gs[4, 4:6])
fig2 = plt.figure(figsize=(14, 15))  
gs2 = gridspec.GridSpec(5, 6, height_ratios=[1, 0.3, 1, 0.05, 1], hspace=0.1, wspace=0.4)   
ax2 = [0, 0, 0, 0, 0,0,0]

ax2[6] = fig2.add_subplot(gs[0, :])  # Top row spanning both columns
ax2[0] = fig2.add_subplot(gs[2, 0:2])  # First subplot in second row (after space)
ax2[1] = fig2.add_subplot(gs[2, 2:4])  # Second subplot in second row
ax2[2] = fig2.add_subplot(gs[2, 4:6])  # First subplot in third row
ax2[3] = fig2.add_subplot(gs[4, 0:2])  # First subplot in third row
ax2[4] = fig2.add_subplot(gs[4, 2:4]) 
ax2[5] = fig2.add_subplot(gs[4, 4:6])
fig.subplots_adjust(left=0.1, right=0.91,top=0.95,bottom=0.08)
fig2.subplots_adjust(left=0.1, right=0.91,top=0.95,bottom=0.08)
ax[6].text(-0.1,0.9,"a)",horizontalalignment='center',verticalalignment='center',transform=ax[6].transAxes,fontsize=28)
ax[0].text(-0.1,0.9,"b)",horizontalalignment='center',verticalalignment='center',transform=ax[0].transAxes,fontsize=28)
ax[3].text(-0.1,0.9,"c)",horizontalalignment='center',verticalalignment='center',transform=ax[3].transAxes,fontsize=28)
ax[1].text(-0.1,0.9,"d)",horizontalalignment='center',verticalalignment='center',transform=ax[1].transAxes,fontsize=28)
ax[4].text(-0.1,0.9,"e)",horizontalalignment='center',verticalalignment='center',transform=ax[4].transAxes,fontsize=28)
ax[2].text(-0.1,0.9,"f)",horizontalalignment='center',verticalalignment='center',transform=ax[2].transAxes,fontsize=28)
ax[5].text(-0.1,0.9,"g)",horizontalalignment='center',verticalalignment='center',transform=ax[5].transAxes,fontsize=28)

ax2[6].text(-0.1,0.9,"a)",horizontalalignment='center',verticalalignment='center',transform=ax2[6].transAxes,fontsize=28)
ax2[0].text(-0.1,0.9,"b)",horizontalalignment='center',verticalalignment='center',transform=ax2[0].transAxes,fontsize=28)
ax2[3].text(-0.1,0.9,"c)",horizontalalignment='center',verticalalignment='center',transform=ax2[3].transAxes,fontsize=28)
ax2[1].text(-0.1,0.9,"d)",horizontalalignment='center',verticalalignment='center',transform=ax2[1].transAxes,fontsize=28)
ax2[4].text(-0.1,0.9,"e)",horizontalalignment='center',verticalalignment='center',transform=ax2[4].transAxes,fontsize=28)
ax2[2].text(-0.1,0.9,"f)",horizontalalignment='center',verticalalignment='center',transform=ax2[2].transAxes,fontsize=28)
ax2[5].text(-0.1,0.9,"g)",horizontalalignment='center',verticalalignment='center',transform=ax2[5].transAxes,fontsize=28)

ax[0].set_title(r"$55-85\degree$")
ax2[0].set_title(r"$55-85\degree$")
ax[2].set_title(r"$110-140\degree$")
ax2[2].set_title(r"$110-140\degree$")
ax[1].set_title(r"$85-110\degree$")
ax2[1].set_title(r"$85-110\degree$")




handles=[]
labels=["Recorded","GLAD-M25","S40RTS"]
colors=[(0.3, 0.3, 0.3),"red","mediumturquoise"]
for i,f in enumerate(["_real","_glad","_8112"]):
    if f=="_real":
        color=(0.285, 0.285, 0.285)
        label="Recorded"
        capsize=8
        elinewidth=2
        fmt="*"
        s=14
        alpha=1
        alpha2=1
    elif f=="_8112":
        color="red"
        label="S40RTS"
        capsize=8
        elinewidth=2
        fmt="."
        s=14
        alpha=1
        alpha2=1
    elif f=="_glad":
        color="dodgerblue"
        label="GLAD-M25"
        capsize=8
        elinewidth=2
        fmt="."
        s=14
        alpha=1
        alpha2=1
   
    
   
    handles.append(Line2D([0], [0], color=color, ls="", marker=fmt, label=labels[i],markersize=12))
                                            
    
    for k in range(0,3):
        
        df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
        df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
            #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
        df_plot=df[(df["0"].isin(ev_list))]
        df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    
        print(len(df_plot["1"]))
        df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
        ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
        ~(df_plot["1"].str.split(".").str[0].isin(us)) &
        ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
        (df_plot["1"].str.split(".").str[0].isin(gsn))]
        df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
        df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<140)]
                # ax[2+k].set_title(r"$100-140\degree$")
                # ax2[2+k].set_title(r"$100-140\degree$")
        for i,df in enumerate(file_names):
                for w in [0,1]:
                    df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
                    df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
                    df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]
        if k==0:  
            #df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
            bin_edges=np.arange(55,146,10)
            bin_centers = bin_edges[:-1] + 5
            
            
            

            df_plot["bin"] = np.digitize(df_plot["7"], bin_edges)
            df_plot["diff1"] = df_plot[str(14+1*5+2)+f] - df_plot[str(14+0*5+2)+f]
            df_plot["diff2"] = df_plot[str(14+1*5+4)+f] - df_plot[str(14+0*5+4)+f]

            # Compute mean and std deviation per bin
            mean_diff1 = df_plot.groupby("bin")["diff1"].mean()
            std_diff1 = df_plot.groupby("bin")["diff1"].std()

            mean_diff2 = df_plot.groupby("bin")["diff2"].mean()
            std_diff2 = df_plot.groupby("bin")["diff2"].std()

            print(len(mean_diff1), bin_centers)

            bin_centers = bin_edges[:-1] + 5
            if f == "_real":
                ax[6].errorbar(bin_centers, mean_diff1, yerr=std_diff1, fmt=fmt, color=color, label=label, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)
                ax2[6].errorbar(bin_centers, mean_diff2, yerr=std_diff2, fmt=fmt, color=color, label=label, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)
            elif f == "_8112":
                ax[6].errorbar(bin_centers-1, mean_diff1, yerr=std_diff1, fmt=fmt, color=color, label=label, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)
                ax2[6].errorbar(bin_centers-1, mean_diff2, yerr=std_diff2, fmt=fmt, color=color, label=label, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)
            elif f == "_glad":
                ax[6].errorbar(bin_centers+1, mean_diff1, yerr=std_diff1, fmt=fmt, color=color, label=label, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)
                ax2[6].errorbar(bin_centers+1, mean_diff2, yerr=std_diff2, fmt=fmt, color=color, label=label, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)
            # else:
            #     ax[4].errorbar(bin_centers, mean_diff1, yerr=std_diff1, fmt=fmt, color=color, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)
            #     ax2[4].errorbar(bin_centers, mean_diff2, yerr=std_diff2, fmt=fmt, color=color, capsize=capsize,capthick=1,elinewidth=elinewidth,ms=s,alpha=alpha)

                
            ax[-1].set_ylabel(r"$\delta \rm T^{ss-s}$")
            ax2[-1].set_ylabel(r"$\delta \rm lnA^{ss/s}$")
            ax[-1].set_xlabel(r"Epicentral distance($\degree$)")
            ax2[-1].set_xlabel(r"Epicentral distance($\degree$)")
            ax[-1].grid("True",which="both",linestyle="dotted")
            ax2[-1].grid("True",which="both",linestyle="dotted")
            ax[-1].set_xlim(50,145)
            ax2[-1].set_xlim(50,145)
            ax[-1].set_ylim(-12,12)
            ax2[-1].set_ylim(-1,1)

            ax_new=ax2[-1].twinx()
            ax_new.set_ylabel(r"$\delta \rm A^{ss/s}$")
            ax_new.set_yticks(ax2[-1].get_yticks())  # Match primary y-ticks
            ax_new.set_yticklabels([exp_format(y, None) for y in ax2[-1].get_yticks()])
        
        if k==0:
            df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<85)]
        elif k==1:
            df_plot=df_plot[(df_plot["7"]>=85) & (df_plot["7"]<110)]
        elif k==2:
            df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<=140)]
           
        
        if 1==1:
            for h in range(2):
                if h==0:
                    col=str(14+1*5+2)
                    axs=ax[0+k]
                    ylabel=r"$\delta \rm T^{ss}$"
                elif h==1:
                    col=str(14+0*5+2)
                    axs=ax[3+k]
                    ylabel=r"$\delta \rm T^{s}$"
                    if k==1:
                       ylabel=r"$\delta \rm T^{s}$"
                if f!="_glad":
                   axs.scatter(df_plot[str(14+1*5+2)+f]-df_plot[str(14+0*5+2)+f],df_plot[col+f],alpha=alpha2,marker=fmt,color=color,label=label,s=8)
                cc=stats.pearsonr(df_plot[str(14+1*5+2)+f]-df_plot[str(14+0*5+2)+f],df_plot[col+f])
                
                if f=="_real":
                   text=f"{round(cc[0], 2):.2f}"
                   axs.text(0.9,0.18,text.ljust(2, '0'),horizontalalignment='center',verticalalignment='center',transform=axs.transAxes,fontsize=20,color=color)
                elif f=="_8112":
                    text=f"{round(cc[0], 2):.2f}"
                    axs.text(0.9,0.05,text.ljust(2, '0'),horizontalalignment='center',verticalalignment='center',transform=axs.transAxes,fontsize=20,color=color)
                elif f=="_glad":
                    text=f"{round(cc[0], 2):.2f}"
                    axs.text(0.9,0.12,text.ljust(2, '0'),horizontalalignment='center',verticalalignment='center',transform=axs.transAxes,fontsize=20,color=color)
                # if k==02
                #    axs.set_xlabel(r"$\delta \rm T^{ss-s}$")
                # elif k==1:
                #     axs.set_xlabel(r"$\delta \rm T^{ss-s}$")
                if h==1:
                    axs.set_xlabel(r"$\delta \rm T^{ss-s}$")
                if k==0:
                    axs.set_ylabel(ylabel)
                axs.set_xlim(-12,12)
                axs.set_ylim(-12,12)
                axs.set_xticks(np.arange(-12,13,4))
                axs.set_yticks(np.arange(-12,13,4))
                axs.yaxis.set_minor_locator(ticker.MultipleLocator(2))
                axs.xaxis.set_minor_locator(ticker.MultipleLocator(2))
                axs.grid("True",which="both",linestyle="dotted")

            for h in range(2):
                if h==0:
                    col=str(14+1*5+4)
                    axs2=ax2[0+k]
                    
                    ylabel=r"$\delta \rm lnA^{ss}$"
                elif h==1:
                    col=str(14+0*5+4)
                    axs2=ax2[3+k]
                    ylabel=r"$\delta \rm lnA^{s}$"
                    if k==1:
                        ylabel=r"$\delta \rm lnA^{s}$"
                if f!="_8112":
                   axs2.scatter(df_plot[str(14+1*5+4)+f]-df_plot[str(14+0*5+4)+f],df_plot[col+f],alpha=alpha2,marker=fmt,color=color,label=label,s=8)
                
                cc=stats.pearsonr(df_plot[str(14+1*5+4)+f]-df_plot[str(14+0*5+4)+f],df_plot[col+f])
                
                if f=="_real":
                   text=f"{round(cc[0], 2):.2f}"
                   axs2.text(0.9,0.18,text.ljust(2, '0'),horizontalalignment='center',verticalalignment='center',transform=axs2.transAxes,fontsize=20,color=color)
                elif f=="_8112":
                    text=f"{round(cc[0], 2):.2f}"
                    axs2.text(0.9,0.05,text.ljust(2, '0'),horizontalalignment='center',verticalalignment='center',transform=axs2.transAxes,fontsize=20,color=color)
                elif f=="_glad":
                    text=f"{round(cc[0], 2):.2f}"
                    axs2.text(0.9,0.12,text.ljust(2, '0'),horizontalalignment='center',verticalalignment='center',transform=axs2.transAxes,fontsize=20,color=color)
                if h==1:
                    axs2.set_xlabel(r"$\delta \rm lnA^{ss/s}$")
                # elif h==1:
                #     axs2.set_xlabel(r"$\delta \rm lnA^{SS/S}$")
                if k==0:    
                   axs2.set_ylabel(ylabel)
                axs2.set_xlim(-2,2)
                axs2.set_ylim(-2,2)
                axs2.set_xticks(np.arange(-2,2.1,1))
                axs2.set_yticks(np.arange(-2,2.1,1))
                axs2.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
                axs2.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
                axs2.grid("True",which="both",linestyle="dotted")
ax[6].legend(handles=handles,loc="lower left",fontsize=20)
ax2[6].legend(handles=handles,loc="lower left",fontsize=20)
tr=np.linspace(-12,12,100)
ar=np.linspace(-2,2,100)

for i in range(6):
    ax[i].plot(tr,tr,"k--")
    ax2[i].plot(ar,-ar,"k--")
    ax2[i].plot(ar,ar,"k--")
    ax[i].plot(tr,-tr,"k--")



#fig.subplots_adjust(top=0.95, bottom=0.1)
fig.savefig(plot_dir+"/relative/scatters/scatter_cc_s40.pdf",dpi=600)
#fig2.subplots_adjust(top=0.95, bottom=0.1)
fig2.savefig(plot_dir+"/relative/scatters/scatter_amp_glad.pdf",dpi=600)
fig.tight_layout()
fig2.tight_layout()
plt.close(fig)
plt.close(fig2)



    
        
fig,ax=plt.subplots(2,2,figsize=(15,15))
gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1])
fig2,ax2=plt.subplots(2,2,figsize=(15,15))
gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 1])
fig_n,ax_n=plt.subplots(1,2,figsize=(15,8))


fig.text(0.02, 0.95, 'a)', fontsize=28, fontweight='bold')  # Position of "a)"
fig.text(0.02, 0.5, 'c)', fontsize=28, fontweight='bold') 
fig.text(0.5, 0.95, 'b)', fontsize=28, fontweight='bold')
fig.text(0.5, 0.5, 'd)', fontsize=28, fontweight='bold')

fig2.text(0.02, 0.92, 'a)', fontsize=28, fontweight='bold')  # Position of "a)"
fig2.text(0.02, 0.5, 'c)', fontsize=28, fontweight='bold') 
fig2.text(0.5, 0.92, 'b)', fontsize=28, fontweight='bold')
fig2.text(0.5, 0.5, 'd)', fontsize=28, fontweight='bold')

fig_n.text(0.02, 0.92, 'a)', fontsize=28, fontweight='bold')  # Position of "a)"
fig_n.text(0.5, 0.92, 'b)', fontsize=28, fontweight='bold')

ax_n[0].set_ylabel(r"$\delta \rm lnA^{ss-s}$-S40")
ax_n[1].set_ylabel(r"$\delta \rm lnA^{ss-s}$-S40")
ax_n[0].set_xlabel(r"$\delta \rm lnA^{ss-s}$-GLAD")
ax_n[1].set_xlabel(r"$\delta \rm lnA^{ss-s}$-GLAD")
for k in range(0,2):
    df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
    df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
        #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))]
    df_plot=df_plot[(df_plot["0"].isin(ev_list))]
    #df_plot=df_plot[~df_plot["1"].str.split(".").str[0].isin(gsn)]
    
    df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    (df_plot["1"].str.split(".").str[0].isin(gsn))]
    df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
    if k==0:
        df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<110)]
        print(len(df_plot["1"]))
        ax[0,k].set_title(r"$ss-s$")
        ax2[0,k].set_title(r"$ss/s$")
    elif k==1:
        df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<=140)]
        print(len(df_plot["1"]))
        ax[1,k].set_title(r"$ss-s$")
        ax2[1,k].set_title(r"$ss/s$")
    elif k==2:
        df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]#df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]
    if k==0:
        ax_n[0].scatter(df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"],df_plot[str(14+1*5+4)+"_8112"]-df_plot[str(14+0*5+4)+"_8112"],alpha=0.5,marker=".",color="red")
        #cc=stats.pearsonr(df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"],df_plot[str(14+1*5+4)+"_8112"]-df_plot[str(14+0*5+4)+"_8112"])
        #ax_n[0].text(0.1,0.9,str(round(cc[0],2)),horizontalalignment='center',verticalalignment='center',transform=ax_n[0].transAxes,fontsize=20,color="red")
    elif k==1:
        ax_n[1].scatter(df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"],df_plot[str(14+1*5+4)+"_8112"]-df_plot[str(14+0*5+4)+"_8112"],alpha=0.5,marker=".",color="red")
        #cc=stats.pearsonr(df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"],df_plot[str(14+1*5+4)+"_8112"]-df_plot[str(14+0*5+4)+"_8112"])
        #ax_n[1].text(0.1,0.9,str(round(cc[0],2)),horizontalalignment='center',verticalalignment='center',transform=ax_n[1].transAxes,fontsize=20,color="red")
    

    for h in range(2):
        if h==0:
            col=str(14+1*5+2)
            axs=ax[0,k]
            ylabel=r"$\delta \rm T^{ss-s}$-S40"
            if k==1:
                ylabel=r"$\delta \rm T^{ss-s}$-S40"
                
        elif h==1:
            col=str(14+0*5+2)
            axs=ax[1,k]
            ylabel=r"$\delta \rm T^{ss-s}$-GLAD"
            if k==1:
                ylabel=r"$\delta \rm T^{ss-s}$-GLAD"
            #ax_n[1].scatter(df_plot[str(14+1*5+2)+"_glad"]-df_plot[str(14+0*5+2)+"_glad"],df_plot[str(14+1*5+2)+"_8112"]-df_plot[str(14+0*5+2)+"_8112"],alpha=0.5,marker=".",color="dodgerblue")
        if h==0:
           axs.scatter(df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"],df_plot[str(14+1*5+2)+"_8112"]-df_plot[str(14+0*5+2)+"_8112"],alpha=0.5,marker=".",color="red")
           cc=stats.pearsonr(df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"],df_plot[str(14+1*5+2)+"_8112"]-df_plot[str(14+0*5+2)+"_8112"])
           axs.text(0.1,0.9,"CC: "+str(round(cc[0],2)),horizontalalignment='center',verticalalignment='center',transform=axs.transAxes,fontsize=20,color="red")

           
        elif h==1:
            axs.scatter(df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"],df_plot[str(14+1*5+2)+"_glad"]-df_plot[str(14+0*5+2)+"_glad"],alpha=0.5,marker=".",color="dodgerblue")
            cc=stats.pearsonr(df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"],df_plot[str(14+1*5+2)+"_glad"]-df_plot[str(14+0*5+2)+"_glad"])
            axs.text(0.1,0.8,"CC: "+str(round(cc[0],2)),horizontalalignment='center',verticalalignment='center',transform=axs.transAxes,fontsize=20,color="dodgerblue")
        if k==0:
            axs.set_xlabel(r"$\delta \rm T^{ss-s}$")
        elif k==1:
            axs.set_xlabel(r"$\delta \rm T^{ss-s}$")
        axs.set_ylabel(ylabel,fontsize=20)
        axs.set_xlim(-12,12)
        axs.set_ylim(-12,12)
        axs.set_xticks(np.arange(-12,13,4))
        axs.set_yticks(np.arange(-12,13,4))
        axs.yaxis.set_minor_locator(ticker.MultipleLocator(2))
        axs.xaxis.set_minor_locator(ticker.MultipleLocator(2))
        axs.grid("True",which="both",linestyle="dotted")

    for h in range(2):
        if h==0:
            col=str(14+1*5+2)
            axs2=ax2[0,k]
            ylabel=r"$\delta \rm lnA^{ss/s}$-S40"
            if k==1:
                ylabel=r"$\delta \rm lnA^{ss/s}$-S40"
        elif h==1:
            col=str(14+0*5+2)
            axs2=ax2[1,k]
            ylabel=r"$\delta \rm lnA^{ss/s}$-GLAD"
            if k==1:
                ylabel=r"$\delta \rm lnA^{ss/s}$-GLAD"
        if h==0:
           axs2.scatter(df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"],df_plot[str(14+1*5+4)+"_8112"]-df_plot[str(14+0*5+4)+"_8112"],alpha=0.5,marker=".",color="red")
           cc=stats.pearsonr(df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"],df_plot[str(14+1*5+4)+"_8112"]-df_plot[str(14+0*5+4)+"_8112"])
           axs2.text(0.1,0.9,str(round(cc[0],2)),horizontalalignment='center',verticalalignment='center',transform=axs2.transAxes,fontsize=20,color="red")
        elif h==1:
            axs2.scatter(df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"],df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"],alpha=0.5,marker=".",color="dodgerblue")
            cc=stats.pearsonr(df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"],df_plot[str(14+1*5+4)+"_glad"]-df_plot[str(14+0*5+4)+"_glad"])
            axs2.text(0.1,0.8,str(round(cc[0],2)),horizontalalignment='center',verticalalignment='center',transform=axs2.transAxes,fontsize=20,color="dodgerblue")
        if k==0:
            axs2.set_xlabel(r"$\delta \rm lnA^{ss/s}$")
        elif k==1:
            axs2.set_xlabel(r"$\delta \rm lnA^{ss/s}$")
        axs2.set_ylabel(ylabel)
        axs2.set_xlim(-2,2)
        axs2.set_ylim(-2,2)
        axs2.set_xticks(np.arange(-2,2.1,0.8))
        axs2.set_yticks(np.arange(-2,2.1,0.8))
        axs2.yaxis.set_minor_locator(ticker.MultipleLocator(0.4))
        axs2.xaxis.set_minor_locator(ticker.MultipleLocator(0.4))
        axs2.grid("True",which="both",linestyle="dotted")

tr=np.linspace(-12,12,100)
ar=np.linspace(-2,2,100)

for i in range(2):
    ax_n[i].plot(ar,ar,"k--")
    ax_n[i].plot(ar,-ar,"k--")
    ax_n[i].set_xlim(-2,2)
    ax_n[i].set_ylim(-2,2)
    for j in range(2):
        ax[i,j].plot(tr,tr,"k--")
        ax2[i,j].plot(ar,-ar,"k--")
        ax2[i,j].plot(ar,ar,"k--")
        ax[i,j].plot(tr,-tr,"k--")
    

fig_n.tight_layout()
fig_n.savefig(plot_dir+"/relative/scatters/scatter_glad_s40.pdf",dpi=600)
#fig.tight_layout()
fig.savefig(plot_dir+"/relative/scatters/"+str(k)+"_cc_real.pdf",dpi=600)
fig2.tight_layout()
fig2.savefig(plot_dir+"/relative/scatters/"+str(k)+"_amp_real.pdf",dpi=600)
plt.close(fig)
plt.close(fig2)  


            
        

  
file_names=[]
labels=[]
n_list=[872,1482,2405,3326,3941,4534,4922,5405,6026,8112,11289,12181,13332,15821,15988]
for n in n_list:
    labels.append(str(n))
    file_names.append("_"+str(n))
file_names.append("_glad")
file_names.append("_real")
file_names.append("_mod_3") 

for k in range(0,2):
    k=ks[k]
    df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
    df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
   #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))] 
    if k==0:
       df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<110)]
    elif k==1:
         df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]
    elif k==2:
         df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]#df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]

   
    df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    (df_plot["1"].str.split(".").str[0].isin(gsn))]
    df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
    print(len(df_plot))
    
    # labels=["PREM","mPREM","S40RTS"]
    # file_names2=["_prem","_mod_3","_8112"]
    # for i,df in enumerate(file_names2):
    #     plt.figure(1,figsize=(15,10))
    #     # if k==0:
    #     plt.scatter(df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df],np.log(df_plot["8"+df]),alpha=0.5,marker="*",color="black")
    #     # else:
    #     #     plt.scatter(df_plot["7"],df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df],alpha=0.5,marker="*",color="black")
        

    #     plt.grid()
    #     plt.ylabel("SS-S Amplitude Ratios")
    #     plt.xlabel("SS-S CC Differnetial Traveltime")
    #     plt.savefig(plot_dir+"/relative/scatters/cc_amp_"+file_names2[i]+"_"+str(k)+"_scatter.png",dpi=600)
    #     plt.close()

    
    pool = mp.Pool(7)
    results=pool.starmap(points, [(rows["3"],rows["2"],rows["4"],rows["6"],rows["5"],rows["7"],2/5.0,3/5.0)  for i,rows in df_plot.iterrows()])
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

    
    
    
    file_names2=["_real","_8112","_glad"]
    labels=["Recorded","S40RTS","GLAD-M25"]
    # labels=["PREM","mPREM","S40RTS"]
    # file_names2=["_prem","_mod_3","_8112"]
    #for ev in ev_list:
        # df_plots=df_plot[df_plot["0"]==ev]  
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    for i,df in enumerate(file_names2):

        
        y,edges=np.histogram(np.log(df_plot["8"+df]),bins=15,range=(-2,2))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            ax.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3,ls="--")
        elif i==1:
            ax.stairs(y/np.sum(y), edges, color="blue",label=labels[i],fill=False,lw=4)
        else:
            ax.stairs(y/np.sum(y), edges, color="red",label=labels[i],fill=False,lw=4)
        
    ax.set_xlabel(r"$ln(R_{SS/S})$")
    ax.set_ylabel("Normalized Frequency")
    if k==2:
       plt.legend(loc="upper right")
    
    
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.025))
    ax.grid(True, which='both', linestyle='dotted', linewidth=1)
    ax.set_xlim(-2,2)
    ax.set_ylim(0,0.2)
    fig.tight_layout()
    plt.savefig(plot_dir+"/relative/histograms/amp_"+str(k)+".pdf",dpi=600)
    plt.close() 

    fig,ax=plt.subplots(1,1,figsize=(10,10))
    for i,df in enumerate(file_names2):

        
        y,edges=np.histogram(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df],bins=15,range=(-2,2))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            ax.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3)
            p=np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df])
        elif i==1:
            ax.stairs(y/np.sum(y), edges, color="blue",label=labels[i],ls="--",fill=False,lw=3)
            mps=np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df])
        else:
            ax.stairs(y/np.sum(y), edges, color="red",label=labels[i],ls="--",fill=False,lw=3)
            mps40=np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df])
    #ax.axvline(p,color="black",ls="--",lw=2)
    #ax.axvline(mps,color="blue",ls="dotted",lw=2) 
    #ax.axvline(mps40,color="red",ls="dotted",lw=1)  
    if k==0:
       ax.set_xlabel(r"$\delta \rm lnA^{ss/s}$")
    else:
        ax.set_xlabel(r"$\delta \rm lnA^{ss/s}$")
    ax.set_ylabel("Normalized Frequency")
    if k==1:
       plt.legend(loc="upper right")
    
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.grid(True, which='both', linestyle='dotted', linewidth=1)
    ax.set_xlim(-2,2)
    fig.tight_layout()
    plt.savefig(plot_dir+"/relative/histograms/dlnA_"+str(k)+".pdf",dpi=600)
    plt.close() 

    fig,ax=plt.subplots(1,1,figsize=(10,10))
    for i,df in enumerate(file_names2):

        
        y,edges=np.histogram(df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df],bins=15,range=(-12,12))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            ax.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3)
            p=np.mean(df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df])
        elif i==1:
            ax.stairs(y/np.sum(y), edges, color="blue",label=labels[i],ls="--",fill=False,lw=3)
            mps=np.mean(df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df])
        else:
            ax.stairs(y/np.sum(y), edges, color="red",label=labels[i],ls="--",fill=False,lw=3)
            mps40=np.mean(df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df])
    ax.axvline(p,color="black",ls="--",lw=2)
    ax.axvline(mps,color="blue",ls="dotted",lw=2) 
    ax.axvline(mps40,color="red",ls="dotted",lw=1)  
    ax.set_xlabel(r"$\delta \rm T^{ss-s}$")
    ax.set_ylabel("Normalized Frequency")
    if k==1:
       plt.legend(loc="upper right")
    
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.grid(True, which='both', linestyle='dotted', linewidth=1)
    ax.set_xlim(-12,12)
    fig.tight_layout()
    plt.savefig(plot_dir+"/relative/histograms/dcc_"+str(k)+".pdf",dpi=600)
    plt.close() 

    
    
    plt.figure(1,figsize=(10,10))
    for i,df in enumerate(file_names2):

        
        y,edges=np.histogram(df_plot[str(14+0*5+4)+df],bins=15,range=(-2,2))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            plt.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3,ls="--")
        elif i==1:
            plt.stairs(y/np.sum(y), edges, color="blue",label=labels[i],fill=False,lw=4)
        else:
            plt.stairs(y/np.sum(y), edges, color="red",label=labels[i],fill=False,lw=4)
        
    plt.xlabel(r"$\delta lnA^s$")
    plt.ylabel("Normalized Frequency")
    if k==1:
       plt.legend(loc="upper right")
    plt.grid()
    plt.tight_layout()
    plt.savefig(plot_dir+"/relative/histograms/amp_S"+str(k)+".pdf",dpi=600)
    plt.close()  

    plt.figure(1,figsize=(10,10))
    for i,df in enumerate(file_names2):

        
        y,edges=np.histogram(df_plot[str(14+1*5+4)+df],bins=15,range=(-2,2))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            plt.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3,ls="--")
        elif i==1:
            plt.stairs(y/np.sum(y), edges, color="blue",label=labels[i],fill=False,lw=4)
        else:
            plt.stairs(y/np.sum(y), edges, color="red",label=labels[i],fill=False,lw=4)
        
    plt.xlabel(r"$\delta lnA^{ss}$")
    plt.ylabel("Normalized Frequency")
    if k==1:
       plt.legend(loc="upper right")
    plt.grid()
    plt.tight_layout()
    plt.savefig(plot_dir+"/relative/histograms/amp_SS"+str(k)+".pdf",dpi=600)
    plt.close()    

    plt.figure(1,figsize=(10,10))
    for i,df in enumerate(file_names2):

        
        y,edges=np.histogram(df_plot[str(14+0*5+2)+df],bins=15,range=(-12,12))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            plt.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3,ls="--")
        elif i==1:
            plt.stairs(y/np.sum(y), edges, color="blue",label=labels[i],fill=False,lw=4)
        else:
            plt.stairs(y/np.sum(y), edges, color="red",label=labels[i],fill=False,lw=4)
        
    plt.xlabel(r"$\delta T^s$")
    plt.ylabel("Normalized Frequency")
    if k==1:
       plt.legend(loc="upper right")
    plt.grid()
    plt.tight_layout()
    plt.savefig(plot_dir+"/relative/histograms/dt_S"+str(k)+".pdf",dpi=600)
    plt.close()  

    plt.figure(1,figsize=(10,10))
    for i,df in enumerate(file_names2):

        
        y,edges=np.histogram(df_plot[str(14+1*5+2)+df],bins=15,range=(-12,12))
        centers = 0.5 * (edges[1:] + edges[:-1])
        if i==0:
            plt.stairs(y/np.sum(y), edges, color="black",label=labels[i],fill=False,lw=3,ls="--")
        elif i==1:
            plt.stairs(y/np.sum(y), edges, color="blue",label=labels[i],fill=False,lw=4)
        else:
            plt.stairs(y/np.sum(y), edges, color="red",label=labels[i],fill=False,lw=4)
        
    plt.xlabel(r"$\delta T^{ss}$")
    plt.ylabel("Normalized Frequency")
    if k==1:
       plt.legend(loc="upper right")
    plt.grid()
    plt.tight_layout()
    plt.savefig(plot_dir+"/relative/histograms/dt_SS"+str(k)+".pdf",dpi=600)
    plt.close()

    print(str(np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df]))+" k="+str(k)+" df")
    for i,df in enumerate(file_names2):
        print(str(np.mean(df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df]))+" k="+str(k)+" "+df)
        fig=pygmt.Figure()
    
        if k==2:
            pygmt.makecpt(cmap='vik',reverse=True,series=[-2,2,0.3],continuous=False) 
        else:
            pygmt.makecpt(cmap='vik',reverse=True,series=[-2,2,0.4],continuous=False)
        df_plot["plot"] = np.log(df_plot["8"+df])
        df_plot["plot"] = df_plot["plot"].apply(lambda n: 2 if n > 2  else -2 if n < -2 else n)
        fig.basemap(region="d",projection="R-150/12c",frame="f")
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.5p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.5p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.2c",fill=df_plot["plot"],pen="0.000000000001p",cmap=True, transparency=40)
        fig.coast(shorelines="1/0.5p",frame="f")
        #if i==0:
        #   fig.colorbar()
        fig.savefig(plot_dir+"/relative/maps/amp_"+str(k)+"_"+str(df)+"_receiver.pdf",dpi=600)
        plt.close()

        fig=pygmt.Figure()
        if k==2:
            pygmt.makecpt(cmap='/home/ayon/WORK/plots_final/col_am.cpt',series=[-1.8,1.8],reverse=True,continuous=False) 
        else:
            pygmt.makecpt(cmap='/home/ayon/WORK/plots_final/col_am.cpt',series=[-1.8,1.8],reverse=True,continuous=False)
        df_plot["plot"] = df_plot[str(14+1*5+4)+df]-df_plot[str(14+0*5+4)+df]
        df_plot["plot"] = df_plot["plot"].apply(lambda n: -1.8 if n > 1.8  else -1.8 if n < -1.8 else n)
        fig.basemap(region="d",projection="R-150/12c",frame="f")
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.4p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.4p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.25c",fill=df_plot["plot"],cmap=True, transparency=40)
        fig.coast(shorelines="1/0.75p",frame="f")
        #if i==0:
        #   fig.colorbar()
        fig.savefig(plot_dir+"/relative/maps/amp_2_"+str(k)+"_"+str(df)+"_receiver_1.5.pdf",dpi=600)
        plt.close()

        fig=pygmt.Figure()
        if k==2:
            pygmt.makecpt(cmap='/home/ayon/WORK/plots_final/col_cc.cpt',series=[-10.5,10.5],continuous=False) 
        else:
            pygmt.makecpt(cmap='/home/ayon/WORK/plots_final/col_cc.cpt',series=[-10.5,10.5],continuous=False)
        df_plot["plot"] = df_plot[str(14+1*5+2)+df]-df_plot[str(14+0*5+2)+df]
        df_plot["plot"] = df_plot["plot"].apply(lambda n: 10 if n > 10  else -10 if n < -10 else n)
        fig.basemap(region="d",projection="R-150/12c",frame="f")
        fig.plot(data=df_plot[["start_lon","start_lat","mid_lon","mid_lat"]],style="=0.01c+s",pen="0.4p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)            
        fig.plot(data=df_plot[["mid_lon","mid_lat","end_lon","end_lat"]],style="=0.01c+s",pen="0.4p,+z",transparency=40)#zvalue=df_plot[16+k*5+4],cmap=True,transparency=40)
        fig.plot(x=df_plot.mid_lon,y=df_plot.mid_lat,style="c0.25c",fill=df_plot["plot"],cmap=True, transparency=40)
        fig.coast(shorelines="1/0.75p",frame="f")
        #if i==0:
        #   fig.colorbar()
        fig.savefig(plot_dir+"/relative/maps/cc_"+str(k)+"_"+str(df)+"_receiver.pdf",dpi=600)
        plt.close()
for k in range(0,2):
    df=pd.read_csv("20_40_SS_S_all.txt",skipinitialspace=True,delimiter=",")
    df_bad=pd.read_csv("./../bad_data_all.txt",skipinitialspace=True,delimiter=",")
        #df_plot=df[df[str(18+k*11+5)+"_real"].notna()]
    df_plot=df[(df["0"].isin(ev_list))]
    if k==0:
        df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<110)]
    #     ax[0,k].set_title(r"$55-100\degree$")
    #     ax2[0,k].set_title(r"$55-100\degree$")
    elif k==1:
         df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<=140)]
    #     ax[0,k].set_title(r"$100-140\degree$")
    #     ax2[0,k].set_title(r"$100-140\degree$")
    # elif k==2:
    #     df_plot=df_plot[(df_plot["7"]>=110) & (df_plot["7"]<140)]#df_plot=df_plot[(np.abs(df_plot[str(18+k*11+2)+"_real"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_1D"])<15) &(np.abs(df_plot[str(18+k*11+2)+"_3D"])<15) & (np.abs(df_plot[str(18+k*11+2)+"_prem_3D"])<15)]
    #df_plot=df_plot[(df_plot["7"]>=55) & (df_plot["7"]<140)]
    df_plot=df_plot[(~((df_plot["5"] >= -135) & (df_plot["5"] <= -59) & (df_plot["6"] >= 32) & (df_plot["6"] <= 80)) &
    ~((df_plot["5"] >= -166) & (df_plot["5"] <= -132) & (df_plot["6"] >= 54) & (df_plot["6"] <= 71)) &
    ~(df_plot["1"].str.split(".").str[0].isin(us)) &
    ~((df_plot["5"] >= -103) & (df_plot["5"] <= -59) & (df_plot["6"] >= 25) & (df_plot["6"] <= 32))) |
    (df_plot["1"].str.split(".").str[0].isin(gsn))]
    df_plot = df_plot[~df_plot[['0', '1']].apply(tuple, axis=1).isin(df_bad[['0', '1']].apply(tuple, axis=1))]
    for i,df in enumerate(file_names):
        for w in [0,1]:
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+2)+df])<12)]
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+1)+df])>0.8)] 
            df_plot=df_plot[(np.abs(df_plot[str(14+w*5+4)+df])<2)]
        
    fig,ax=plt.subplots(1,1,figsize=(8,8))
    #fig.text(0.02, 0.92, 'a)', fontsize=28, fontweight='bold')  # Position of "a)"
    #fig.text(0.02, 0.5, 'c)', fontsize=28, fontweight='bold') 
    #fig.text(0.5, 0.92, 'b)', fontsize=28, fontweight='bold')
    #fig.text(0.5, 0.5, 'd)', fontsize=28, fontweight='bold')

    #for i,f in enumerate(file_names2):
    #ax[0].scatter(df_plot["7"],df_plot[str(14+1*5+2)+"_real"]-df_plot[str(14+0*5+2)+"_real"],alpha=0.5,marker=".",color="blue")
    #ax[0].scatter(df_plot["7"],df_plot[str(14+1*5+2)+"_mod_3"]-df_plot[str(14+0*5+2)+"_mod_3"],alpha=0.5,marker=".",color="red")
    ax.scatter(df_plot["7"],df_plot[str(14+1*5+4)+"_real"]-df_plot[str(14+0*5+4)+"_real"],alpha=0.5,marker=".",color="blue")
    #ax[1].scatter(df_plot["7"],df_plot[str(14+1*5+4)+"_mod_3_2"]-df_plot[str(14+0*5+4)+"_mod_3_2"],alpha=0.5,marker=".",color="red")
    ax.set_xlabel(r"Epicentral Distaces$\degree$")

    #ax[0].set_ylabel(r"$\delta \rm T^{SS-S}$")
    ax.set_ylabel(r"$\delta \rm lnA^{ss-s}$")
    #ax[0].set_xlim(55,140)
    #ax[1].set_xlim(55,140)
    #ax[0].set_ylim(-12,12)
    ax.set_ylim(-2,2)
    #ax[0].set_xticks(np.arange(55,141,10))
    #ax[1].set_xticks(np.arange(55,141,10))
    #ax[0].set_yticks(np.arange(-12,13,4))
    ax.set_yticks(np.arange(-2,2.1,0.5))
    #ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.2))

    fig.tight_layout()
    fig.savefig(plot_dir+"/relative/scatters/SS_S_distance"+str(k)+".pdf",dpi=600)
    plt.close(fig)

        

    # df_np=df_plot[(df_plot["mid_lon"] < -105 ) & (df_plot["mid_lat"] <0) & (df_plot["mid_lat"] > -45)] 
    # df_np.to_csv("regions/"+str(k)+"SS_S_south_pacific.txt",index=False,na_rep=np.nan)
    # df_si=df_plot[(df_plot["mid_lon"] > 90 ) & (df_plot["mid_lon"] <115) & (df_plot["mid_lat"] < -30) & (df_plot["mid_lat"] > -60)]   
    # df_si.to_csv("regions/"+str(k)+"SS_S_south_indian.txt",index=False,na_rep=np.nan)
    # df_np=df_plot[(df_plot["mid_lon"] <-115 ) & (df_plot["mid_lon"] > -135) & (df_plot["mid_lat"] > 60) & (df_plot["mid_lat"] < 90)]  
    # df_np.to_csv("regions/"+str(k)+"SS_S_arctic.txt",index=False,na_rep=np.nan)
    # df_sa=df_plot[(df_plot["mid_lon"] <0 ) & (df_plot["mid_lon"] >-45) & (df_plot["mid_lat"] <0) & (df_plot["mid_lat"] > -60)]
    # df_sa.to_csv("regions/"+str(k)+"SS_S_south_atlantic.txt",index=False,na_rep=np.nan)
    # df_npc=df_plot[(df_plot["mid_lon"] > -175 ) & (df_plot["mid_lon"] <-120) & (df_plot["mid_lat"] > 0) & (df_plot["mid_lat"] < 45)]
    # df_npc.to_csv("regions/"+str(k)+"SS_S_north_pacific.txt",index=False,na_rep=np.nan)
    # df_al=df_plot[((df_plot["mid_lon"] > 150 ) | (df_plot["mid_lon"] < -155)) & (df_plot["mid_lat"] > 45) & (df_plot["mid_lat"] < 65)]
    # df_al.to_csv("regions/"+str(k)+"SS_S_alaska.txt",index=False,na_rep=np.nan)
    # df_np=df_plot[(df_plot["mid_lon"] > 105 ) & (df_plot["mid_lon"] < 120) & (df_plot["mid_lat"] > -45) & (df_plot["mid_lat"] < -15)]
    # df_np.to_csv("regions/"+str(k)+"SS_S_west_aus.txt",index=False,na_rep=np.nan)
    # df_np=df_plot[(df_plot["mid_lon"] < -15 ) & (df_plot["mid_lon"] > -45) & (df_plot["mid_lat"] > 45) & (df_plot["mid_lat"] < 75)]
    # df_np.to_csv("regions/"+str(k)+"SS_S_north_atlantic.txt",index=False,na_rep=np.nan)
    # df_np=df_plot[(df_plot["mid_lon"] > -105 ) & (df_plot["mid_lon"] < -90) & (df_plot["mid_lat"] > -15) & (df_plot["mid_lat"] < 15)]
    # df_np.to_csv("regions/"+str(k)+"SS_S_parallel_south.txt",index=False,na_rep=np.nan)
    # df_np=df_plot[(df_plot["mid_lon"] > -120 ) & (df_plot["mid_lon"] < -105) & (df_plot["mid_lat"] > 0) & (df_plot["mid_lat"] < 30)]
    # df_np.to_csv("regions/"+str(k)+"SS_S_parallel_north.txt",index=False,na_rep=np.nan)

        
            
    
    

