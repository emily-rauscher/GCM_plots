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

def igcm_olr(path, runname, oom, surfp, radea, createplot, savefig, savenamelw,savenamesw, BOTH,LastOrb):
    # wave == 'LW' (1) or "SW' (2)
    # both == total energy outgoing (SW + LW)
    
    if LastOrb==False:     
        ##### load files #####
        with open(path+runname+'/fort.64') as f:
            first_line=f.readline()
            nlat,nlon,nlev=first_line.split()
            nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)
            print '  '
            print ' ....reading fort.64 (LW)'
            print '       nlat=', nlat, 'nlon=', nlon, 'nlev=', nlev

        #if nlat!=len(lat_arr) or nlon!=len(lon_arr):
        #    print 'ERROR: Lat or Lon mis-match, check files'
        #    exit()

        data_64=np.empty([nlat*nlon,3])*0.0
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
        l=0
        while l<data_64.shape[0]:
            lon_arr=np.append(lon_arr,lon_arr_f[l])
            l+=nlat

        lat_arr=data_64[:nlat,1]

        data_lw=np.empty([nlon,nlat,3])
        for l in range(0,data_64.shape[0]):
            lon,lat=data_64[l,:2]
            lon_i,lat_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0]
            data_lw[lon_i,lat_i,:]=data_64[l,:]
            
        data_lw=data_lw[:,:,2]

        print 'LW', data_lw.shape

        #######################
        if BOTH==True:
            with open(path+runname+'/fort.65') as f:
                first_line=f.readline()
                nlat,nlon,nlev=first_line.split()
                nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)
                print '  '
                print ' ....reading fort.65 (SW)'
                print '       nlat=', nlat, 'nlon=', nlon, 'nlev=', nlev

            #if nlat!=len(lat_arr) or nlon!=len(lon_arr):
            #    print 'ERROR: Lat or Lon mis-match, check files'
            #    exit()

            data_65=np.empty([nlat*nlon,3])*0.0
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
            l=0
            while l<data_65.shape[0]:
                lon_arr=np.append(lon_arr,lon_arr_f[l])
                l+=nlat

            lat_arr=data_65[:nlat,1]

            data_sw=np.empty([nlon,nlat,3])
            for l in range(0,data_65.shape[0]):
                lon,lat=data_65[l,:2]
                lon_i,lat_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0]
                data_sw[lon_i,lat_i,:]=data_65[l,:]
            
            data_sw=data_sw[:,:,2]
            print 'SW', data_sw.shape
            
    if LastOrb==True:
        print ' DOING LAST ORBIT AVERAGES....'
        filen=6400
        while filen<6490:
            file_lw='fort.'+str(int(filen))
            if filen==6400:
                with open(path+runname+'/'+file_lw) as f:
                    first_line=f.readline()
                    nlat,nlon,nlev=first_line.split()
                    nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)
                    data_64=np.empty([nlat*nlon])*0.0
                    print '  '
                    print ' ....reading fort.6400 (LW)'
                    print '       nlat=', nlat, 'nlon=', nlon, 'nlev=', nlev
                    print '  '
            if filen>6400:
                print ' ....reading fort.',filen,' (LW)'

            l=0
            data_64h=np.empty([nlat*nlon,3])*0.0
            with open(path+runname+'/'+file_lw) as f:
                for line in f:
                    if l==0:
                        l+=1
                        continue
                    if l>nlat*nlon:
                        continue
                        l+=1
                    #print line, line.split()
                    data_64h[l-1] = line.split()
                    l+=1
                print '       END OF FILE: DONE'
            f.close()
            
            lon_arr_f=data_64h[:,0]
            lon_arr=np.array([])
            l=0
            while l<data_64h.shape[0]:
                lon_arr=np.append(lon_arr,lon_arr_f[l])
                l+=nlat

            lat_arr=data_64h[:nlat,1]

            data_lwh=np.empty([nlon,nlat,3])

            for l in range(0,data_64h.shape[0]):
                lon,lat=data_64h[l,:2]
                #print l,lon,lat
                lon_i,lat_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0]
                #print l,lon,lat,lon_i,lat_i

                data_lwh[lon_i,lat_i,:]=data_64h[l,:]
            
            data_lwh=data_lwh[:,:,2]
            
            if filen==6400:
                data_lw=np.copy(data_lwh)
            elif filen>6400:
                data_lw=np.dstack((data_lw,data_lwh))
            filen+=1
                             

            #print '            (', data_64.shape, ' )'
        print data_lw.shape
        data_lw=np.nanmean(data_lw,axis=2)
        print data_lw.shape

        print 'LW', data_lw.shape

        #######################
        if BOTH==True:
            filen=6500
            while filen<6590:
                file_sw='fort.'+str(int(filen))
                if filen==6500:
                    with open(path+runname+'/'+file_sw) as f:
                        first_line=f.readline()
                        nlat,nlon,nlev=first_line.split()
                        nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)
                        data_65=np.empty([nlat*nlon,3])*0.0
                        print '  '
                        print ' ....reading fort.6500 (SW)'
                        print '       nlat=', nlat, 'nlon=', nlon, 'nlev=', nlon
                        print '  '
                if filen>6500:
                    print  ' ....reading fort.',filen,' (SW)'

                l=0
                data_65h=np.empty([nlat*nlon,3])*0.0
                with open(path+runname+'/'+file_sw) as f:
                    for line in f:
                        if l==0:
                            l+=1
                            continue
                        if l>nlat*nlon:
                            continue
                            l+=1
                        #print line, line.split()
                        data_65h[l-1] = line.split()
                        l+=1
                    print '       END OF FILE: DONE'
                f.close()
                
                lon_arr_f=data_65h[:,0]
                lon_arr=np.array([])
                l=0
                while l<data_65h.shape[0]:
                    lon_arr=np.append(lon_arr,lon_arr_f[l])
                    l+=nlat

                lat_arr=data_65h[:nlat,1]

                data_swh=np.empty([nlon,nlat,3])
                for l in range(0,data_65h.shape[0]):
                    lon,lat=data_65h[l,:2]
                    lon_i,lat_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0]
                    data_swh[lon_i,lat_i,:]=data_65h[l,:]
                
                data_swh=data_swh[:,:,2]
               
                if filen==6500:
                    data_sw=np.copy(data_swh)
                elif filen>6500:
                    data_sw=np.dstack((data_sw,data_swh))
                filen+=1

            print data_sw.shape
            data_sw=np.nanmean(data_sw,axis=2)
            print data_sw.shape
            
            print 'SW', data_sw.shape

    
    #######################
    # Sum values for total outgoing radiation (taken from igcm_olr.pro)
    # ' Total integrated output (W): ',total(olr*cos(lat*!pi/180.))*(2.*!pi/nlon)*(!pi/nlat)*radea^2
    total_lw=np.nansum(data_lw[:,:]*np.cos(lat_arr*np.pi/180.))*(2*np.pi/nlon)*(np.pi/nlat)*radea**2.
    dayside_lw=((np.nansum(data_lw[0:np.int(nlon/4),:]*np.cos(lat_arr*np.pi/180.))  
                                      +np.nansum(data_lw[np.int(nlon*3./4):nlon-1,:]*np.cos(lat_arr*np.pi/180.))) 
                                    *(2.*np.pi/nlon)*(np.pi/nlat)*radea**2)
    if BOTH==True:
        total_sw=np.nansum(data_sw[:,:]*np.cos(lat_arr*np.pi/180.))*(2*np.pi/nlon)*(np.pi/nlat)*radea**2.
        dayside_sw=((np.nansum(data_sw[0:np.int(nlon/4),:]*np.cos(lat_arr*np.pi/180.)) 
                                      +np.nansum(data_sw[np.int(nlon*3./4):nlon-1,:]*np.cos(lat_arr*np.pi/180.))) 
                                    *(2.*np.pi/nlon)*(np.pi/nlat)*radea**2)
    print '******************************'
    print 'Total Integrated Output (W):'
    print '  LW:', total_lw
    if BOTH==True:
        print '  SW:', total_sw
        print ' sum:', total_lw+total_sw
    print '-------------------------------'
    print ' Dayside Integrated Output (W):'
    print '  LW:', dayside_lw
    if BOTH==True:
        print '  SW:', dayside_sw
        print ' sum:', dayside_lw+dayside_sw
    print '******************************'
    
    #######################
    if createplot==True:
        
        plt.figure(1,figsize=(10,6.25))
        #shifting arrays so centered at lat and lon=0.0
        CENTER=np.zeros([nlon,nlat])
        for i in range(0,len(lon_arr)):
            if i<nlon/2:
                CENTER[i+nlon/2,:]=data_lw[i,:]
            if i>=nlon/2:
                CENTER[i-nlon/2,:]=data_lw[i,:]
        plt_data=CENTER
        plt_lon=np.linspace(-180,180,nlon)
        
        LON,LAT=np.meshgrid(plt_lon,lat_arr)
        
        cbar_levs=np.round_(np.linspace(np.nanmin(plt_data)/1.01,np.nanmax(plt_data)*1.01,20),2)
        
        p=plt.contourf(LON,LAT,plt_data.T,levels=cbar_levs,cmap=plt.cm.Reds,zorder=0)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18)
        
        plt.ylim(np.nanmin(lat_arr),np.nanmax(lat_arr))
        plt.xlim(np.nanmin(plt_lon),np.nanmax(plt_lon))
        
        plt.ylabel('Latitude [${\circ}$]',fontsize=20)
        plt.xlabel('Longitude [${\circ}$]',fontsize=20)
        plt.yticks(fontsize=18,fontproperties=font)
        plt.xticks(fontsize=18,fontproperties=font)
        
        plt.figtext(0.77,0.5,'LW Outgoing Radiation (W)',
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
                    CENTER[i+nlon/2,:]=data_sw[i,:]
                if i>=nlon/2:
                    CENTER[i-nlon/2,:]=data_sw[i,:]
            plt_data=CENTER
            plt_lon=np.linspace(-180,180,nlon)

            LON,LAT=np.meshgrid(plt_lon,lat_arr)

            cbar_levs=np.round_(np.linspace(np.nanmin(plt_data)/1.01,np.nanmax(plt_data)*1.01,20),2)
            if np.all(np.diff(cbar_levs) > 0):
                p=plt.contourf(LON,LAT,plt_data.T,levels=cbar_levs,cmap=plt.cm.Blues,zorder=0)
            else:
                print '   ...calc levels not working, autogen...'
                p=plt.contourf(LON,LAT,plt_data.T,cmap=plt.cm.Blues,zorder=0)
            c=plt.colorbar(p)
            c.ax.tick_params(labelsize=18)

            plt.ylim(np.nanmin(lat_arr),np.nanmax(lat_arr))
            plt.xlim(np.nanmin(plt_lon),np.nanmax(plt_lon))

            plt.ylabel('Latitude [${\circ}$]',fontsize=20)
            plt.xlabel('Longitude [${\circ}$]',fontsize=20)
            plt.yticks(fontsize=18,fontproperties=font)
            plt.xticks(fontsize=18,fontproperties=font)

            plt.figtext(0.77,0.5,'SW Outgoing Radiation (W)',
                        fontsize=20,rotation='vertical',ha='center',va='center')

            if savefig==True:
                plt.savefig(savenamesw,rasterized=True,transparent=True)
            plt.show()
    if BOTH==True:
        return data_lw, data_sw, total_lw, total_sw,lon_arr,lat_arr
    else:
        return data_lw,np.nan, total_lw,np.nan,lon_arr,lat_arr