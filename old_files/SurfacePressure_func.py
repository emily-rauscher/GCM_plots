from load_data_noin import load_data
import numpy as np

import pickle

import matplotlib
import matplotlib.pyplot as plt
#%matplotlib inline
from matplotlib.font_manager import FontProperties
font0=FontProperties()
font=font0.copy()
font.set_family('serif')

fontb=font.copy()
fontb.set_weight('bold')

params = {'font.family': 'serif',}
matplotlib.rcParams.update(params)

def SurfPress(data_50, oom, surfp, createplot, savefig, savename):
    # load fort.26
    #runname,lon_arr,lat_arr,oom,surfp,p_BAR,data_26,data_lo,data_olr=load_data(path,runname,oom,surfp,False,False,verbose)

    # load fort.50
    
    if createplot==True:
        
        plt.figure(figsize=(10,6.25))
        #shifting arrays so centered at lat and lon=0.0
        CENTER=np.zeros([nlon,nlat])
        for i in range(0,len(lon_arr)):
            if i<nlon/2:
                CENTER[i+nlon/2,:]=data_50[i,:,2]
            if i>=nlon/2:
                CENTER[i-nlon/2,:]=data_50[i,:,2]
        plt_data=CENTER/surfp
        plt_lon=np.linspace(-180,180,nlon)
        
        LON,LAT=np.meshgrid(plt_lon,lat_arr)
        
        cbar_levs=np.round_(np.linspace(np.nanmin(plt_data)/1.01,np.nanmax(plt_data)*1.01,20),2)
        
        if np.all(np.diff(cbar_levs) > 0):
            p=plt.contourf(LON,LAT,plt_data.T,levels=cbar_levs,cmap=plt.cm.Purples,zorder=0)
        else:
            p=plt.contourf(LON,LAT,plt_data.T,cmap=plt.cm.Purples,zorder=0)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18)
        
        plt.ylim(np.nanmin(lat_arr),np.nanmax(lat_arr))
        plt.xlim(np.nanmin(plt_lon),np.nanmax(plt_lon))
        
        plt.ylabel('Latitude [${\circ}$]',fontsize=20)
        plt.xlabel('Longitude [${\circ}$]',fontsize=20)
        plt.yticks(fontsize=18,fontproperties=font)
        plt.xticks(fontsize=18,fontproperties=font)
        
        plt.figtext(0.77,0.5,'Relative Surface Pressure',fontsize=20,rotation='vertical',ha='center',va='center')
    
        if savefig==True:
            plt.savefig(savename,rasterized=True,transparent=True)
        plt.show()
    return data_50