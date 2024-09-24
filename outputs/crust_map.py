import pygmt
import numpy as np
import matplotlib.pyplot as plt

evlos=[179.29]
evlas=[-23.29]
stlos=[13.0011]
stlas=[46.4142]
plot_dir="/scratch1/09038/ayon8181/synthetic/scripts_amp/outputs/today"
data=np.loadtxt('/scratch1/09038/ayon8181/scripts_amp/outputs/depthtomoho.xyz')
data[:, 2] = data[:, 2] * -1
#crust=pygmt.xyz2grd(data,region="d",spacing=(0.5,0.5))

fig=pygmt.Figure()
fig.basemap(region="d",projection="N-70/12c",frame=True)
grid = pygmt.surface(data, region="d", spacing=[0.5,0.5], verbose=True)
fig.grdimage(grid=grid,cmap='polar',projection="N-70/12c",verbose="q")
fig.colorbar(frame="x+lCrustal Thickness(km)")
for i in range(len(evlos)):
    stlo=stlos[i]
    stla=stlas[i]
    evlo=evlos[i]
    evla=evlas[i]
    fig.plot(x=[stlo,evlo],y=[stla,evla],pen="0.5p")
    fig.plot(x=stlo,y=stla,style="i0.2c",fill='green')
    fig.plot(x=stlo,y=stla,style="c0.2c",fill='yellow')
fig.coast(shorelines=True, frame=True)
        #if i==0:
        #   fig.colorbar()
fig.savefig(plot_dir+"/map.png",dpi=600)
plt.close()