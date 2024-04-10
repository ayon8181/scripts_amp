#!/usr/bin/env python
import pyasdf
import obspy
import numpy as np
from scipy import integrate
from scipy.signal import hilbert
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy.signal.cross_correlation import correlate_template
from obspy.imaging import beachball
from scipy import fft
from pyasdf import ASDFDataSet
import matplotlib.pyplot as plt
import time
import sys
import math
import csv
"""
event       = str(sys.argv[1])
obsd_tag    = str(sys.argv[2])
obsd_fname  = str(sys.argv[3])
"""

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
    az=[]
    np_1=[]
    for i in range(1,360):
        rp=np.abs(rad_pat_sh(np.deg2rad(rake),np.deg2rad(dip),np.deg2rad(i),np.deg2rad(strike),np.deg2rad(20.0)))
        if rp<=0.001:
           np_0.append(i)       
    return np_0

def calulate_SNR


ds=ASDFDataSet("/work2/09038/ayon8181/frontera/new_events/synthetics/asdf_data/202206042338A.real_data.h5",mode="r")
event=ds.events[0]

print(nodal_plane_filter(event))


# Read Events

# Calculate Nodal Planes for SH

# Remove if Azimuth is in +/- 15 degree

# Calculate SNR for Energy

# Calculate SNR for Power Spectra

#