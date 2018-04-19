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

def Temp_Plot(lons,lats,press,data,lev,latcenter,loncenter,temp_ex,units_a,units_t,savefig,savename,ortho,ver,cbarL,cbarM):
    lon_arr=lons
    lat_arr=lats
    nlon=len(lon_arr)
    nlat=len(lat_arr)
    
    p_BAR=press
    
    data_26=data
    #######################################
    colors1 = plt.cm.YlGnBu_r(np.linspace(0, 1, 128))
    colors2 = plt.cm.YlOrBr(np.linspace(0., 1, 128))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    #########################################
    if units_a==0:
        lat_arr=lat_arr*(np.pi/180.)
        lon_arr=lon_arr*(np.pi/180.)
    if units_a==1:
        lat_arr=lat_arr
        lon_arr=lon_arr
        
    if units_t==1:
        data_26-=273.15
    if units_t==2:
        data_26=(9./5.)*(data_26-273.15)+32.
        
    minT=np.int((np.nanmin(data_26)))*1.0
    maxT=np.ceil(np.nanmax(data_26))
    
    if cbarL>0 and cbarM>0:
        if units_t==1:
            cbarL-=273.15
            cbarM-=273.15
        if units_t==2:
            cbarL=(9./5.)*(cbarL-273.15)+32.
            cbarM=(9./5.)*(cbarM-273.15)+32.
        minT=cbarL
        maxT=cbarM
    cbar_levs=np.linspace(minT-temp_ex,maxT+temp_ex,(maxT-minT+2*temp_ex)+1)
    
    if ver==True:
        print '-------------------------------------------------------'
        print '  Plotting Atmosphere Level: ', lev, '-> Pressure=',np.round(p_BAR[lev],4)
        print '-------------------------------------------------------'
        if units_t==1:
            print 'Min Temp [C], Plot limit: ', np.nanmin(data_26), minT
            print 'Max Temp [C], Plot limit: ', np.nanmax(data_26), maxT
        if units_t==2:
            print 'Min Temp [F], Plot limit: ', np.nanmin(data_26), minT
            print 'Max Temp [F], Plot limit: ', np.nanmax(data_26), maxT
        else:
            print 'Min Temp [K], Plot limit: ', np.nanmin(data_26), minT
            print 'Max Temp [K], Plot limit: ', np.nanmax(data_26), maxT
        print '-------------------------------------------------------'
    
    if ortho==False:
    
        plt.figure(figsize=(10,6.25))
        plt.gcf().subplots_adjust(bottom=0.12,top=0.95,left=0.12,right=0.99)

        #shifting arrays so centered at lat and lon=0.0
        TEMPS=np.zeros([nlon,nlat])
        for i in range(0,len(lon_arr)):
            if i<nlon/2:
                TEMPS[i+nlon/2,:]=data_26[i,:]
            if i>=nlon/2:
                TEMPS[i-nlon/2,:]=data_26[i,:]
        plt_data=TEMPS
        plt_lon=np.linspace(-180,180,nlon)
      

        LON,LAT=np.meshgrid(plt_lon,lat_arr)
        
        p=plt.contourf(LON,LAT,plt_data.T,levels=cbar_levs,cmap=mymap,zorder=0)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18) 

        plt.grid(color='white',linewidth=0.5,linestyle='--',alpha=0.5, zorder=10)

        if units_t==0:
            plt.figtext(0.84,0.5,'Temperature [K]',fontsize=20,rotation='vertical',ha='center',va='center')
        if units_t==1:
            plt.figtext(0.84,0.5,'Temperature [$^{\circ}$ C]',fontsize=20,rotation='vertical',ha='center',va='center')
        if units_t==2:
            plt.figtext(0.84,0.5,'Temperature [$^{\circ}$ F]',fontsize=20,rotation='vertical',ha='center',va='center')

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
        
    if ortho==True:
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

        fig=plt.figure(figsize=(11,8))
        ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        ax.axis('off')
        
        p=plt.contourf(X,Y,plt_data0*MASK,levels=cbar_levs,cmap=mymap,zorder=0)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18) 
        
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
        if units_t==0:
            plt.figtext(0.8,0.5,'Temperature [K]',fontsize=20,rotation='vertical',ha='center',va='center')
        if units_t==1:
            plt.figtext(0.8,0.5,'Temperature [$^{\circ}$ C]',fontsize=20,rotation='vertical',ha='center',va='center')
        if units_t==2:
            plt.figtext(0.8,0.5,'Temperature [$^{\circ}$ F]',fontsize=20,rotation='vertical',ha='center',va='center')

        #plt.yticks(fontsize=18,fontproperties=font)
        #plt.xticks(fontsize=18,fontproperties=font)

      
        if savefig==True:
            plt.savefig(savename,rasterized=True)
        if ver==True:
            plt.show()
        
        plt.close()
        
    return
