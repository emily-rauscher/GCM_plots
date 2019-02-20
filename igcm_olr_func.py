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

def igcm_olr(path, runname, oom, surfp, createplot, savefig, savenamesw,savenaelw, WAVE, BOTH):
    # wave == 'LW' (1) or "SW' (2)
    # both == total energy outgoing (SW + LW)
    
    if WAVE==1:
        file='fort.64'
    if wave==2:
        file='fort.65'
        
    ##### load files #####
    with open(path+runname+'/fort.64') as f:
        first_line=f.readline()
        nlat,nlon=first_line.split()
        nlat,nlon=int(nlat),int(nlon)
        print '  '
        print ' ....reading fort.64 (LW)'
        print '       nlat=', nlat, 'nlon=', nlon
    
    #if nlat!=len(lat_arr) or nlon!=len(lon_arr):
    #    print 'ERROR: Lat or Lon mis-match, check files'
    #    exit()
        
    data64=np.empty([nlat*nlon,3])*0.0
    l=0
    with open(path+runname+'/fort.64') as f:
        for line in f:
            if l==0:
                l+=1
                continue
            if l>nlat*nlon:
                continue
                l+=1
            #print line, line.split()
            data_64[l-1] = line.split()
            l+=1
        print '       END OF FILE: DONE'
    f.close()
    
    print data_64.shape
    
    lon_arr_f=data_64[:,0]
    lon_arr=np.array([])
    for l in range(0,len(lon_arr_f)):
        el=lon_arr_f[l]
        if not el in lon_arr:
            lon_arr=np.append(lon_arr,el)

    lat_arr_f=data_64[:,1]
    lat_arr=np.array([])
    for l in range(0,len(lat_arr_f)):
        el=lat_arr_f[l]
        if not el in lat_arr:
            lat_arr=np.append(lat_arr,el)
            
    data_lw=np.empty([nlon,nlat,3])
    for l in range(0,data_64.shape[0]):
        lon,lat=data_64[l,:2]
        lon_i,lat_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0]
        data_lw[lon_i,lat_i,:]=data_64[l,:]
            
    print 'LW', data_lw.shape
    
    #######################
    if BOTH==True:
        with open(path+runname+'/fort.65') as f:
            first_line=f.readline()
            nlat,nlon=first_line.split()
            nlat,nlon=int(nlat),int(nlon)
            print '  '
            print ' ....reading fort.65 (SW)'
            print '       nlat=', nlat, 'nlon=', nlon

        #if nlat!=len(lat_arr) or nlon!=len(lon_arr):
        #    print 'ERROR: Lat or Lon mis-match, check files'
        #    exit()

        data65=np.empty([nlat*nlon,3])*0.0
        l=0
        with open(path+runname+'/fort.65') as f:
            for line in f:
                if l==0:
                    l+=1
                    continue
                if l>nlat*nlon:
                    continue
                    l+=1
                #print line, line.split()
                data_65[l-1] = line.split()
                l+=1
            print '       END OF FILE: DONE'
        f.close()

        print data_65.shape

        lon_arr_f=data_65[:,0]
        lon_arr=np.array([])
        for l in range(0,len(lon_arr_f)):
            el=lon_arr_f[l]
            if not el in lon_arr:
                lon_arr=np.append(lon_arr,el)

        lat_arr_f=data_65[:,1]
        lat_arr=np.array([])
        for l in range(0,len(lat_arr_f)):
            el=lat_arr_f[l]
            if not el in lat_arr:
                lat_arr=np.append(lat_arr,el)

        data_sw=np.empty([nlon,nlat,3])
        for l in range(0,data_65.shape[0]):
            lon,lat=data_65[l,:2]
            lon_i,lat_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0]
            data_sw[lon_i,lat_i,:]=data_65[l,:]

        print 'SW', data_sw.shape
    
    #######################
    # Sum values for total outgoing radiation (taken from igcm_olr.pro)
    
    #######################
    if createplot==True:
        
        plt.figure(1,figsize=(10,6.25))
        #shifting arrays so centered at lat and lon=0.0
        CENTER=np.zeros([nlon,nlat])
        for i in range(0,len(lon_arr)):
            if i<nlon/2:
                CENTER[i+nlon/2,:]=data_lw[i,:,2]
            if i>=nlon/2:
                CENTER[i-nlon/2,:]=data_lw[i,:,2]
        plt_data=CENTER
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
        
        plt.figtext(0.77,0.5,'LW Outgoing Radiation',
                    fontsize=20,rotation='vertical',ha='center',va='center')
    
        if savefig==True:
            plt.savefig(savenamelw,rasterized=True,transparent=True)
        plt.show()
        
        if BOTH==True:
            plt.figure(2,figsize=(10,6.25))
            #shifting arrays so centered at lat and lon=0.0
            CENTER=np.zeros([nlon,nlat])
            for i in range(0,len(lon_arr)):
                if i<nlon/2:
                    CENTER[i+nlon/2,:]=data_sw[i,:,2]
                if i>=nlon/2:
                    CENTER[i-nlon/2,:]=data_sw[i,:,2]
            plt_data=CENTER
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

            plt.figtext(0.77,0.5,'SW Outgoing Radiation',
                        fontsize=20,rotation='vertical',ha='center',va='center')

            if savefig==True:
                plt.savefig(savenamesw,rasterized=True,transparent=True)
            plt.show()
    if BOTH==True:
        return data_lw, data_sw, total_lw, total_Sw
    else:
        return data_lw, total_lw