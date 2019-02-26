from load_data_noin import load_data
import numpy as np

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

def SurfPress(path, runname, oom, surfp, createplot, savefig, savename):
    # load fort.26
    #runname,lon_arr,lat_arr,oom,surfp,p_BAR,data_26,data_lo,data_olr=load_data(path,runname,oom,surfp,False,False,verbose)

    # load fort.50
    with open(path+runname+'/fort.50') as f:
        first_line=f.readline()
        nlat,nlon=first_line.split()
        nlat,nlon=int(nlat),int(nlon)
        print '  '
        print ' ....reading fort.50'
        print '       nlat=', nlat, 'nlon=', nlon
    
    #if nlat!=len(lat_arr) or nlon!=len(lon_arr):
    #    print 'ERROR: Lat or Lon mis-match, check files'
    #    exit()
        
    data50=np.empty([nlat*nlon,3])*0.0
    l=0
    with open(path+runname+'/fort.50') as f:
        for line in f:
            if l==0:
                l+=1
                continue
            if l>nlat*nlon:
                continue
                l+=1
            #print line, line.split()
            data50[l-1] = line.split()
            l+=1
        print '       END OF FILE: DONE'
    f.close()
    
    print data50.shape
    
    lon_arr_f=data_50[:,0]
    lon_arr=np.array([])
    l=0
    while l<data_50.shape[0]:
        lon_arr=np.append(lon_arr,lon_arr_f[l])
        l+=nlat

    lat_arr=data_50[:nlat,1]
            
    data_50=np.empty([nlon,nlat,3])
    for l in range(0,data50.shape[0]):
        lon,lat=data50[l,:2]
        lon_i,lat_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0]
        data_50[lon_i,lat_i,:]=data50[l,:]
            
    print data_50.shape
    data_50[:,:,2]=(data_50[:,:,2]+1.)*surfp
    
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
        
        p=plt.contourf(LON,LAT,plt_data.T,levels=cbar_levs,cmap=plt.cm.Purples,zorder=0)
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