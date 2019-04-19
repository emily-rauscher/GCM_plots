import numpy as np
import math

from scipy.interpolate import griddata

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches

from matplotlib.font_manager import FontProperties
font0=FontProperties()
font=font0.copy()
font.set_family('serif')

fontb=font.copy()
fontb.set_weight('bold')

params = {'font.family': 'serif',}
matplotlib.rcParams.update(params)

def press_Plot(plot,lons,lats,press,data,units_a,units_t,units_w,freeze,caption,cap,
          savefig,savename,zeros,ver,cbarL,cbarM,cbar_even,ex,ncolors,longavg,long_pl,lo,oom):
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
    
    #data_26 index: Nlev,nlon,nlat,nparam
    
    if plot==0: #TEMPERATURE
        if lo==True:
            ind=0
        else:
            ind=5 
        
        #colors1 = plt.cm.YlGnBu_r(np.linspace(0, 1, 128))
        #colors2 = plt.cm.YlOrBr(np.linspace(0., 1, 128))
        cm_data = np.loadtxt("ScientificColourMaps5/roma/roma.txt")
        colors=np.flip(cm_data,axis=0)
        
    if plot==1:
        if lo==True:
            ind=1
        else:
            ind=3 # uwind
        
        cm_data = np.loadtxt("ScientificColourMaps5/cork/cork.txt")
        cm_data=np.flip(cm_data,axis=0)
        #colors1 = plt.cm.BuGn_r(np.linspace(0, 1, 128))
        #colors2 = plt.cm.RdPu(np.linspace(0., 1, 128))
        colors = cm_data#np.vstack((colors1, colors2))
    if plot==2:
        if lo==True:
            ind=2
        else:
            ind=4 # vwind
        
        colors1 = plt.cm.Purples_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.Oranges(np.linspace(0., 1, 128))
        colors = np.vstack((colors1, colors2))
    
    data_26=np.copy(data[:,:,:,ind])
    #######################################

    # combine them and build a new colormap
    
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
    ##########################################################
    if ver==True:
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
        ##########################################################
    if longavg==True:
        plt_data=np.nanmedian(data_26,axis=1)
    else:
        lon_indx=np.argmin(np.abs(lon_arr-long_pl))
        plt_data=data_26[:,lon_indx,:]

    print plt_data.shape
    plt.figure(figsize=(8.25,8))
    plt.gcf().subplots_adjust(bottom=0.12,top=0.95,left=0.12,right=0.99)
    
    LAT,PRESS_P=np.meshgrid(lat_arr,p_BAR)
    
    p=plt.contourf(LAT,PRESS_P,plt_data,levels=cbar_levs,cmap=mymap,zorder=0,alpha=1.0)
    c=plt.colorbar(p)
    c.ax.tick_params(labelsize=18)
    if freeze==True:
        plt.contour(PRESS_P,LAT,plt_data,levels=[freezeT],colors='black',zorder=1)
    if zeros==True:
        zerolevels=0.0
        plt.contour(LAT,PRESS_P,plt_data,levels=[zerolevels],colors='black',zorder=1)

    plt.grid(color='white',linewidth=0.5,linestyle='--',alpha=0.5, zorder=10)
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

    plt.ylim(np.nanmax(p_BAR),np.nanmin(p_BAR))
    plt.xlim(np.nanmin(lat_arr),np.nanmax(lat_arr))

    if oom >0:
        plt.yscale('log')
    if units_a==1:
        plt.xlabel('Latitude [degrees]',fontsize=20)
    else:
        plt.xlabel('Latitude [radians]',fontsize=20)

    plt.ylabel('Pressure [bar]',fontsize=20)
    #adding a square caption
    if caption==True:
        rectangle=patches.Rectangle((45,.00001), 205,.000025, fill=True, fc='#ECE0DD')
        plt.gca().add_patch(rectangle)
        plt.annotate(cap, (46,.00002), fontsize='13', color='black', weight='bold')
    if savefig==True:
        plt.savefig(savename,rasterized=True)
    if ver==True:
        plt.show()

    plt.close()
    return

    
