import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from matplotlib.font_manager import FontProperties
font0=FontProperties()
font=font0.copy()
font.set_family('serif')

fontb=font.copy()
fontb.set_weight('bold')

params = {'font.family': 'serif',}
matplotlib.rcParams.update(params)

def ortholine(R,lonl,latl,loncenter,latcenter):
    x=R*np.cos(latl)*np.sin(lonl-loncenter)
    y=R*(np.cos(latcenter)*np.sin(latl)-np.sin(latcenter)*np.cos(latl)*np.cos(lonl-loncenter))
    cosc=np.sin(latcenter)*np.sin(latl)+np.cos(latcenter)*np.cos(latl)*np.cos(lonl-loncenter)
    for i in range(0,len(cosc)):
        if cosc[i]<0:
            x[i]=np.nan
            y[i]=np.nan
    return x,y

def igcm_Plot(plot,lons,lats,press,data,lev,latcenter,loncenter,ex,units_a,units_t,units_w,
              savefig,savename,ortho,ver,cbarL,cbarM,cbar_even,ncolors,vfrac):
    lon_arr=lons
    lat_arr=lats
    nlon=len(lon_arr)
    nlat=len(lat_arr)
    
    p_BAR=press
    
    # nparam index: 
    #      0=lons
    #      1=lats
    #      2=levs
    #      3=u wind
    #      4=v wind
    #      5=temps
    
    if plot==0: #TEMPERATURE
        ind=5 
        
        colors1 = plt.cm.YlGnBu_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.YlOrBr(np.linspace(0., 1, 128))
        
    if plot==1:
        ind=3 # uwind
        
        colors1 = plt.cm.BuGn_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.RdPu(np.linspace(0., 1, 128))
    if plot==2:
        ind=4 # vwind
        
        colors1 = plt.cm.Purples_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.Oranges(np.linspace(0., 1, 128))
        
    
    if plot==3:  #PLOTTING STREAM LINES with faded temps behind
        ind=5  #(want contours to be temp)
        colors1 = plt.cm.YlGnBu_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.YlOrBr(np.linspace(0., 1, 128))
    
        data_26_u=np.copy(data[lev,:,:,3])
        data_26_v=np.copy(data[lev,:,:,4])
        #data_26_W=(np.sqrt(np.square(data_26_u)+np.square(data_26_v)))
        
    data_26=np.copy(data[lev,:,:,ind])
       
    #######################################

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    #########################################
    # set up lat and lon units
    if units_a==0: #radians
        lat_arr=lat_arr*(np.pi/180.)
        lon_arr=lon_arr*(np.pi/180.)
    if units_a==1: #degrees (as input)
        lat_arr=lat_arr
        lon_arr=lon_arr
    
    if plot==0 or plot==3: #TEMPERATURE, convert units    
        if units_t==1:
            data_26-=273.15
        if units_t==2:
            data_26=(9./5.)*(data_26-273.15)+32.
    if plot==1 or plot==2: #WINDS, convert units
        if units_w==1:
            data_26/=1000.  #m/s to km/s
        if units_w==2:
            data_26=(data_26/1609.34)*3600.  #m/s to mph
        
    minV=np.int((np.nanmin(data_26)))*1.0
    maxV=np.ceil(np.nanmax(data_26))
    
    if plot==1 or plot==2:
        if cbar_even==True:
            lim=np.nanmax([np.abs(minV),np.abs(maxV)])
            minV=-1.0*lim
            maxV=lim
    
    if cbarL>0 and cbarM>0:
        if plot==0 or plot==3:
            if units_t==1:
                cbarL-=273.15
                cbarM-=273.15
            if units_t==2:
                cbarL=(9./5.)*(cbarL-273.15)+32.
                cbarM=(9./5.)*(cbarM-273.15)+32.
        if plot==1 or plot==2:
            if units_w==1:
                cbarL/=1000.  #m/s to km/s
                cbarM/=1000.  #m/s to km/s
            if units_w==2:
                cbarL=(cbarL/1609.34)*3600.  #m/s to mph
                cbarM=(cbarM/1609.34)*3600.  #m/s to mph
        minV=cbarL
        maxV=cbarM
    
    if ncolors==0:
        cbar_levs=np.linspace(minV-ex,maxV+ex,(maxV-minV+2*ex)+1)
    else:
        cbar_levs=np.linspace(minV-ex,maxV+ex,ncolors*((maxV-minV+2*ex))+1.0)
    
    if ver==True:
        print '-------------------------------------------------------'
        print '  Plotting Atmosphere Level: ', lev, '-> Pressure=',np.round(p_BAR[lev],4)
        print '-------------------------------------------------------'
        if plot==0:
            if units_t==1:
                print 'Min Temp [C], Plot limit: ', np.nanmin(data_26), minV
                print 'Max Temp [C], Plot limit: ', np.nanmax(data_26), maxV
            if units_t==2:
                print 'Min Temp [F], Plot limit: ', np.nanmin(data_26), minV
                print 'Max Temp [F], Plot limit: ', np.nanmax(data_26), maxV
            if units_t==0:
                print 'Min Temp [K], Plot limit: ', np.nanmin(data_26), minV
                print 'Max Temp [K], Plot limit: ', np.nanmax(data_26), maxV
            print '-------------------------------------------------------'
        if plot==1:
            if units_w==0:
                print 'Min UWind [m/s], Plot limit: ', np.nanmin(data_26), minV
                print 'Max UWind [m/s], Plot limit: ', np.nanmax(data_26), maxV
            if units_w==1:
                print 'Min UWind [km/s], Plot limit: ', np.nanmin(data_26), minV
                print 'Max UWind [km/s], Plot limit: ', np.nanmax(data_26), maxV
            if units_w==2:
                print 'Min UWind [mph], Plot limit: ', np.nanmin(data_26), minV
                print 'Max UWind [mph], Plot limit: ', np.nanmax(data_26), maxV
            print '-------------------------------------------------------'
        if plot==2:
            if units_w==0:
                print 'Min VWind [m/s], Plot limit: ', np.nanmin(data_26), minV
                print 'Max VWind [m/s], Plot limit: ', np.nanmax(data_26), maxV
            if units_w==1:
                print 'Min VWind [km/s], Plot limit: ', np.nanmin(data_26), minV
                print 'Max VWind [km/s], Plot limit: ', np.nanmax(data_26), maxV
            if units_w==2:
                print 'Min VWind [mph], Plot limit: ', np.nanmin(data_26), minV
                print 'Max VWind [mph], Plot limit: ', np.nanmax(data_26), maxV
            print '-------------------------------------------------------'
    
    if ortho==False:
    
        plt.figure(figsize=(10,6.25))
        plt.gcf().subplots_adjust(bottom=0.12,top=0.95,left=0.12,right=0.99)

        #shifting arrays so centered at lat and lon=0.0
        CENTER=np.zeros([nlon,nlat])
        for i in range(0,len(lon_arr)):
            if i<nlon/2:
                CENTER[i+nlon/2,:]=data_26[i,:]
            if i>=nlon/2:
                CENTER[i-nlon/2,:]=data_26[i,:]
        plt_data=CENTER
        plt_lon=np.linspace(-180,180,nlon)
        
        if plot==3:
            CENTER_U=np.zeros([nlon,nlat])
            CENTER_V=np.zeros([nlon,nlat])
            for i in range(0,len(lon_arr)):
                if i<nlon/2:
                    CENTER_U[i+nlon/2,:]=data_26_u[i,:]
                    CENTER_V[i+nlon/2,:]=data_26_v[i,:]
                if i>=nlon/2:
                    CENTER_U[i-nlon/2,:]=data_26_u[i,:]
                    CENTER_V[i-nlon/2,:]=data_26_v[i,:]
            plt_U0=CENTER_U
            plt_V0=CENTER_V
      

        LON,LAT=np.meshgrid(plt_lon,lat_arr)
        
        if plot==3:
            a=0.3
        else:
            a=1.0
        p=plt.contourf(LON,LAT,plt_data.T,levels=cbar_levs,cmap=mymap,zorder=0,alpha=a)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18) 

        plt.grid(color='white',linewidth=0.5,linestyle='--',alpha=0.5, zorder=10)
        
        if plot==3:
            data_26_W=(np.sqrt(np.square(plt_U0)+np.square(plt_V0)))
            lw=5.*data_26_W/np.nanmax(data_26_W)
            plt.streamplot(LON,LAT,plt_U0.T,plt_V0.T,density=vfrac,color=data_26_W.T,cmap=matplotlib.cm.bone_r,linewidth=lw.T)
        
        if plot==0:
            if units_t==0:
                plt.figtext(0.84,0.5,'Temperature [K]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_t==1:
                plt.figtext(0.84,0.5,'Temperature [$^{\circ}$ C]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_t==2:
                plt.figtext(0.84,0.5,'Temperature [$^{\circ}$ F]',fontsize=20,rotation='vertical',ha='center',va='center')
        if plot==1:
            if units_w==0:
                plt.figtext(0.84,0.5,'E-W Wind [m/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==1:
                plt.figtext(0.84,0.5,'E-W Wind [km/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==2:
                plt.figtext(0.84,0.5,'E-W Wind [mph]',fontsize=20,rotation='vertical',ha='center',va='center')
        if plot==2:
            if units_w==0:
                plt.figtext(0.84,0.5,'N-S Wind [m/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==1:
                plt.figtext(0.84,0.5,'N-S Wind [km/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==2:
                plt.figtext(0.84,0.5,'N-S Wind [mph]',fontsize=20,rotation='vertical',ha='center',va='center')

        plt.yticks(fontsize=18,fontproperties=font)
        plt.xticks(fontsize=18,fontproperties=font)

        if units_a==1:
            plt.ylabel('Latitude [${\circ}$]',fontsize=20)
            plt.xlabel('Longitude [${\circ}$]',fontsize=20)
        else:
            plt.ylabel('Latitude [radians]',fontsize=20)
            plt.xlabel('Longitude [radians]',fontsize=20)
      
        if savefig==True:
            plt.savefig(savename,rasterized=True)
        if ver==True:
            plt.show()
        
        plt.close()
        return
        
    elif ortho==True:
        ol=1 
        R=1.0
        if units_a==1:  #convert to rads for ortho projection
            lat_arr=lat_arr*(np.pi/180.)
            lon_arr=lon_arr*(np.pi/180.)
            
        latcenter=latcenter*np.pi/180.
        loncenter=loncenter*np.pi/180.
            
        lon_arr=np.append(lon_arr,lon_arr[:ol])
        #lat_arr=np.append(lat_arr,lat_arr[:ol])
        LON,LAT=np.meshgrid(lon_arr,lat_arr)

        X=np.zeros_like(LON)
        Y=np.zeros_like(LAT)

        MASK=np.ones_like(X)
        

        for i in range(0,X.shape[0]):
            for j in range(0,X.shape[1]):
                X[i,j]=R*np.cos(LAT[i,j])*np.sin(LON[i,j]-loncenter)
                Y[i,j]=R*(np.cos(latcenter)*np.sin(LAT[i,j])-np.sin(latcenter)*np.cos(LAT[i,j])*np.cos(LON[i,j]-loncenter))
                cosc=np.sin(latcenter)*np.sin(LAT[i,j])+np.cos(latcenter)*np.cos(LAT[i,j])*np.cos(LON[i,j]-loncenter)
                if cosc<0:
                    MASK[i,j]=np.nan

        plt_data0=np.zeros_like(X)
        for i in range(0,X.shape[0]):
            for j in range(0,X.shape[1]):
                if i>=nlat and j<nlon:
                    plt_data0[i,j]=(data_26.T)[nlat+ol-i,j]
                if j>=nlon and i <nlat:
                    plt_data0[i,j]=(data_26.T)[i,nlon+ol-j] 
                if j>=nlon and i>=nlat:
                    plt_data0[i,j]=(data_26.T)[nlat+ol-i,nlon+ol-j] 
        plt_data0[:nlat,:nlon]=(data_26.T)
        
        if plot==3:
            plt_U0=np.zeros_like(X)
            plt_V0=np.zeros_like(X)
            for i in range(0,X.shape[0]):
                for j in range(0,X.shape[1]):
                    if i>=nlat and j<nlon:
                        plt_U0[i,j]=(data_26_u.T)[nlat+ol-i,j]
                        plt_V0[i,j]=(data_26_v.T)[nlat+ol-i,j]
                    if j>=nlon and i <nlat:
                        plt_U0[i,j]=(data_26_u.T)[i,nlon+ol-j] 
                        plt_V0[i,j]=(data_26_v.T)[i,nlon+ol-j]
                    if j>=nlon and i>=nlat:
                        plt_U0[i,j]=(data_26_u.T)[nlat+ol-i,nlon+ol-j] 
                        plt_V0[i,j]=(data_26_v.T)[nlat+ol-i,nlon+ol-j]
            plt_U0[:nlat,:nlon]=(data_26_u.T)
            plt_V0[:nlat,:nlon]=(data_26_v.T)
            

        fig=plt.figure(figsize=(11,8))
        ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        ax.axis('off')
        
        if plot==3:
            a=0.3
        else:
            a=1.0
        p=plt.contourf(X,Y,plt_data0*MASK,levels=cbar_levs,cmap=mymap,zorder=0,alpha=a)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18) 
        
        if plot==3:
            data_26_W=(np.sqrt(np.square(plt_U0)+np.square(plt_V0)))
            lw=5.*data_26_W/np.nanmax(data_26_W)
            print X.shape, Y.shape
            print plt_U0.shape, plt_V0.shape
            print data_26_W.shape, lw.shape
            plt.streamplot(X,Y,plt_U0,plt_V0,density=vfrac,color=data_26_W,cmap=matplotlib.cm.bone_r,linewidth=lw)
        
        #lat lines
        x,y=ortholine(R,lon_arr,0*np.pi/180.,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        x,y=ortholine(R,lon_arr,30.*np.pi/180.,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        x,y=ortholine(R,lon_arr,-30.*np.pi/180.,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        x,y=ortholine(R,lon_arr,45.*np.pi/180.,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        x,y=ortholine(R,lon_arr,-45.*np.pi/180.,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        x,y=ortholine(R,lon_arr,60.*np.pi/180.,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        x,y=ortholine(R,lon_arr,-60.*np.pi/180.,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        #lon lines
        x,y=ortholine(R,0*np.pi/180.,lat_arr,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        x,y=ortholine(R,30.*np.pi/180.,lat_arr,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        x,y=ortholine(R,-30.*np.pi/180.,lat_arr,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        x,y=ortholine(R,60.*np.pi/180.,lat_arr,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        x,y=ortholine(R,-60.*np.pi/180.,lat_arr,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        x,y=ortholine(R,90.*np.pi/180.,lat_arr,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        x,y=ortholine(R,-90.*np.pi/180.,lat_arr,loncenter,latcenter)
        plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
        
        
        #colorbar label
        if plot==0:
            if units_t==0:
                plt.figtext(0.8,0.5,'Temperature [K]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_t==1:
                plt.figtext(0.8,0.5,'Temperature [$^{\circ}$ C]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_t==2:
                plt.figtext(0.8,0.5,'Temperature [$^{\circ}$ F]',fontsize=20,rotation='vertical',ha='center',va='center')
        if plot==1:
            if units_w==0:
                plt.figtext(0.8,0.5,'E-W Wind [m/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==1:
                plt.figtext(0.8,0.5,'E-W Wind [km/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==2:
                plt.figtext(0.8,0.5,'E-W Wind [mph]',fontsize=20,rotation='vertical',ha='center',va='center')
        if plot==2:
            if units_w==0:
                plt.figtext(0.8,0.5,'N-S Wind [m/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==1:
                plt.figtext(0.8,0.5,'N-S Wind [km/s]',fontsize=20,rotation='vertical',ha='center',va='center')
            if units_w==2:
                plt.figtext(0.8,0.5,'N-S Wind [mph]',fontsize=20,rotation='vertical',ha='center',va='center')

        #plt.yticks(fontsize=18,fontproperties=font)
        #plt.xticks(fontsize=18,fontproperties=font)

      
        if savefig==True:
            plt.savefig(savename,rasterized=True)
        if ver==True:
            plt.show()
        
        plt.close()
        
        return
