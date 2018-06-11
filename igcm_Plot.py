import numpy as np
import math

from scipy.interpolate import griddata

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

def igcm_Plot(plot,lons,lats,press,data,lev,latcenter,loncenter,ex,units_a,units_t,units_w,freeze,
              savefig,savename,ortho,ver,cbarL,cbarM,cbar_even,ncolors,vfrac,lo):
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
    
    #lo index:
    # 0 =temp
    # 1 = u wind
    # 2 = v wind
    
    
    
    if plot==0: #TEMPERATURE
        if lo==True:
            ind=0
        else:
            ind=5 
        
        colors1 = plt.cm.YlGnBu_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.YlOrBr(np.linspace(0., 1, 128))
        
    if plot==1:
        if lo==True:
            ind=1
        else:
            ind=3 # uwind
        
        colors1 = plt.cm.BuGn_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.RdPu(np.linspace(0., 1, 128))
    if plot==2:
        if lo==True:
            ind=2
        else:
            ind=4 # vwind
        
        colors1 = plt.cm.Purples_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.Oranges(np.linspace(0., 1, 128))
        
    
    if plot==3:  #PLOTTING STREAM LINES with faded temps behind
        if lo==True:
            ind=0
        else:
            ind=5  #(want contours to be temp)
        colors1 = plt.cm.YlGnBu_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.YlOrBr(np.linspace(0., 1, 128))
        
        greys=plt.cm.gray_r(np.linspace(0., 0.75, 128))
        mygreys=mcolors.LinearSegmentedColormap.from_list('my_colormap', greys)
        if lo==True:
            data_26_u=np.copy(data[lev,:,:,1])
            data_26_v=np.copy(data[lev,:,:,2])
        else:
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
    
    freezeT=273.15
    if plot==0 or plot==3: #TEMPERATURE, convert units    
        if units_t==1:
            data_26-=273.15
            freezeT=0.0
        if units_t==2:
            data_26=(9./5.)*(data_26-273.15)+32.
            freezeT=32.0
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
    
    if cbarL!=0 and cbarM!=0:
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
        minV=np.int((np.nanmin(cbarL)))*1.0
        maxV=np.ceil(np.nanmax(cbarM))
    
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
    ###################################
    # CYLINDRICAL PROJECTION PLOTTING #
    ###################################
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
            #a=0.3
            a=0.7 #add a transparency for temps behind streamplot overlay
        else:
            a=1.0
        p=plt.contourf(LON,LAT,plt_data.T,levels=cbar_levs,cmap=mymap,zorder=0,alpha=a)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18)
        if freeze==True:
            plt.contour(LON,LAT,plt_data.T,levels=[freezeT],colors='black',zorder=1)

        plt.grid(color='white',linewidth=0.5,linestyle='--',alpha=0.5, zorder=10)
        
        if plot==3:
            data_26_W=(np.sqrt(np.square(plt_U0)+np.square(plt_V0)))
            lw=5.*data_26_W/np.nanmax(data_26_W)
            plt.streamplot(LON,LAT,plt_U0.T,plt_V0.T,density=vfrac,color='black',linewidth=lw.T,zorder=10)
        
        if plot==0 or plot==3:
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
        
        plt.ylim(np.nanmin(lat_arr),np.nanmax(lat_arr))
        plt.xlim(np.nanmin(plt_lon),np.nanmax(plt_lon))

        if units_a==1:
            plt.ylabel('Latitude [${\circ}$]',fontsize=20)
            plt.xlabel('Longitude [${\circ}$]',fontsize=20)
        else:
            plt.ylabel('Latitude [radians]',fontsize=20)
            plt.xlabel('Longitude [radians]',fontsize=20)
      
        if savefig==True:
            plt.savefig(savename,rasterized=True,transparent=True)
        if ver==True:
            plt.show()
        
        plt.close()
        return
    
    ####################################
    # ORTHOGRAPHIC PROJECTION PLOTTING #
    ####################################
    elif ortho==True:
        ol=1 
        R=1.0
        if units_a==1:  #convert to rads for ortho projection (needed to go to x/y)
            lat_arr=lat_arr*(np.pi/180.)
            lon_arr=lon_arr*(np.pi/180.)
            
        latcenter=latcenter*np.pi/180.
        loncenter=loncenter*np.pi/180.
            
        lon_arr=np.append(lon_arr,lon_arr[:ol]) #extend longitudes to overlap 
        LON,LAT=np.meshgrid(lon_arr,lat_arr)
        
        # Extends plot data to overlap in lon space 
        plt_data0=np.zeros_like(LON)
        for i in range(0,LON.shape[0]):
            for j in range(0,LON.shape[1]):
                if i>=nlat and j<nlon:
                    plt_data0[i,j]=(data_26.T)[nlat+ol-i,j]
                if j>=nlon and i <nlat:
                    plt_data0[i,j]=(data_26.T)[i,nlon+ol-j] 
                if j>=nlon and i>=nlat:
                    plt_data0[i,j]=(data_26.T)[nlat+ol-i,nlon+ol-j] 
        plt_data0[:nlat,:nlon]=(data_26.T)
        
        if plot==3:  #if plotting streamplot, overlap U and V
            plt_U0=np.zeros_like(LON)
            plt_V0=np.zeros_like(LON)
            for i in range(0,LON.shape[0]):
                for j in range(0,LON.shape[1]):
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

        #X=np.zeros_like(LON) # X and Y dummy arrays
        #Y=np.zeros_like(LAT)

        MASK=np.ones_like(LON) # to mask points on unseen hemisphere
        
        X=R*np.cos(LAT)*np.sin(LON-loncenter)
        Y=R*(np.cos(latcenter)*np.sin(LAT)-np.sin(latcenter)*np.cos(LAT)*np.cos(LON-loncenter))
        cosc=np.sin(latcenter)*np.sin(LAT)+np.cos(latcenter)*np.cos(LAT)*np.cos(LON-loncenter)
        
        for i in range(0,X.shape[0]):  #checks for unseen hemisphere
            for j in range(0,X.shape[1]):
                if cosc[i,j]<0:
                    MASK[i,j]=np.nan
  
        #if plot==3:
        #    plt.figure(3)
        #    plt.contourf(LON,LAT,plt_U0)
        #    plt.show()
        #plt.figure(4)
        #plt.contourf(LON,LAT,plt_data0)
        #plt.show()
            
        fig=plt.figure(1,figsize=(11,8))
        ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        ax.axis('off')
        
        if plot==3:
            a=0.7  #add a transparency to the temps for streamplot overlay
        else:
            a=1.0
        p=plt.contourf(X,Y,plt_data0*MASK,levels=cbar_levs,cmap=mymap,zorder=0,alpha=a)
        c=plt.colorbar(p)
        c.ax.tick_params(labelsize=18) 
        
        if freeze==True:
            plt.contour(LON,LAT,plt_data.T,levels=freeze,colors='black',zorder=1)
        
        if plot==3:
            # have to intperolate onto grid, need to flatten first
            
            U0_DATA_flat=(plt_U0*MASK).flatten()
            V0_DATA_flat=(plt_V0*MASK).flatten()
            #temp_flat=(plt_data0*MASK).flatten()
            
            X_nonan_U=(X.flatten())[~np.isnan(U0_DATA_flat)]
            Y_nonan_U=(Y.flatten())[~np.isnan(U0_DATA_flat)]
            u_nonan=(U0_DATA_flat)[~np.isnan(U0_DATA_flat)]
            
            X_nonan_V=(X.flatten())[~np.isnan(V0_DATA_flat)]
            Y_nonan_V=(Y.flatten())[~np.isnan(V0_DATA_flat)]
            v_nonan=(V0_DATA_flat)[~np.isnan(V0_DATA_flat)]
            
            
            new_x=np.linspace(np.nanmin(X),np.nanmax(X),X.shape[0])
            new_y=np.linspace(np.nanmin(Y),np.nanmax(Y),Y.shape[1])
            X1,Y1=np.meshgrid(new_x,new_y)
            # send this to lat lon for interpolation (otherwise ignores back hemisphere)
            rho=np.sqrt(X1*X1 +Y1*Y1)
            c=np.arcsin(rho/R)
            new_lat=np.arcsin(np.cos(c)*np.sin(latcenter)+Y1*np.sin(c)*np.cos(latcenter)/rho)
            new_lon=loncenter+np.arctan2((X1*np.sin(c)),(rho*np.cos(c)*np.cos(latcenter)-Y1*np.sin(c)*np.sin(latcenter)))
            
            interp_U0_cy=griddata((LON.flatten(),LAT.flatten()),plt_U0.flatten(),(new_lon.flatten(),new_lat.flatten()),
                                  fill_value='nan')
            interp_V0_cy=griddata((LON.flatten(),LAT.flatten()),plt_V0.flatten(),(new_lon.flatten(),new_lat.flatten()),
                                  fill_value='nan')
            
            cosc_new=np.sin(latcenter)*np.sin(new_lat)+np.cos(latcenter)*np.cos(new_lat)*np.cos(new_lon-latcenter)
            mask_new=np.ones_like(X1)
            for i in range(0,X1.shape[0]):
                for j in range(0,X1.shape[1]):
                    if cosc_new[i,j]<0.0:
                        mask_new[i,j]=np.nan
            ##############################################
            
            mask1=np.ones_like(X1)
            for i in range(0,X1.shape[0]):
                for j in range(0,X1.shape[1]):
                    r=np.sqrt(X1[i,j]*X1[i,j]+Y1[i,j]*Y1[i,j])
                    if r>1:
                        mask1[i,j]=np.nan
            
            #interp_U0=griddata((X_nonan_U,Y_nonan_U),u_nonan,(X1.flatten(),Y1.flatten()),
            #                   fill_value='nan')
            #interp_V0=griddata((X_nonan_V,Y_nonan_V),v_nonan,(X1.flatten(),Y1.flatten()),
            #                   fill_value='nan')

            #plt_U0_int=interp_U0.reshape(X1.shape)
            #plt_V0_int=interp_V0.reshape(X1.shape)
            
            plt_U0_int=interp_U0_cy.reshape(X1.shape)
            plt_V0_int=interp_V0_cy.reshape(X1.shape)
            
            data_26_W=(np.sqrt(np.square(plt_U0_int)+np.square(plt_V0_int)))
            lw=5.*data_26_W/np.nanmax(data_26_W)
            
            plt.streamplot(X1,Y1,plt_U0_int,plt_V0_int,density=vfrac,color='black',linewidth=lw,zorder=10)
            #plt.contour(X1,Y1,plt_t0_int,levels=cbar_levs,colors='black')
        
        #lat lines for orthoprojection
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
        if plot==0 or plot==3:
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
            plt.savefig(savename,rasterized=True,transparent=True)
        if ver==True:
            plt.show()
        
        plt.close()
        
        #if plot==3:
        #    plt.figure(2)
        #    plt.contourf(new_lon,new_lat,plt_U0_int)
        #    plt.show()
            
        #plt.figure(2)
        #plt.scatter(X,Y,c=plt_U0*MASK)
        #plt.figure(3)
        #plt.scatter(X1,Y1,c=plt_U0_int*mask1)
        
        #plt.figure(4)
        #plt.scatter(X,Y,c=plt_V0*MASK)
        #plt.figure(5)
        #plt.scatter(X1,Y1,c=plt_V0_int*mask1)
        
        return
