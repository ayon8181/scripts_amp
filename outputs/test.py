from obspy.taup import TauPyModel
import numpy as np

# Define a TauP model. You can use existing models like 'prem' or 'ak135'.
model = TauPyModel(model="prem")

rp=model.get_ray_paths_geo(source_depth_in_km=200, source_latitude_in_deg=32, source_longitude_in_deg=23,receiver_latitude_in_deg=-78,receiver_longitude_in_deg=-92,phase_list=['SSS'], resample=True)
print(rp[0].path)