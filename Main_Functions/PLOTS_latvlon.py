import numpy as np
import os

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
#%matplotlib inline
from matplotlib.font_manager import FontProperties

from matplotlib import font_manager as fm, rcParams
params = {'font.family': 'serif',}
matplotlib.rcParams.update(params)

def typical_colorbars(plottype):
    if 'lw' in plottype:
        cm = plt.cm.OrRd
    if 'sw' in plottype:
        cm = plt.cm.PuBu
    if 'uw' in plottype:
        cm = np.loadtxt("ScientificColourMaps5/cork/cork.txt")
        cm = np.flip(cm, axis=0)
    if 'vw' in plottype:
        colors1 = plt.cm.Purples_r(np.linspace(0, 1, 128))
        colors2 = plt.cm.Oranges(np.linspace(0., 1, 128))
        cm      = np.vstack((colors1, colors2))
    if 'tp' in plottype:
        cm = np.loadtxt("ScientificColourMaps5/lajolla/lajolla.txt")
        cm = np.flip(cm, axis=0)
    return cm

def gen_levels(zfl_load,ind,filen,cbar_even,cbarL,cbarM,ncolors,ex):
    minV = np.int((np.nanmin(zfl_load)))*1.0
    maxV = np.ceil(np.nanmax(zfl_load))

    if ind == 3 or ind == 4 and 'fort26.txt' in filen:
        if cbar_even == True:
            lim = np.nanmax([np.abs(minV),np.abs(maxV)])
            minV = -1.0 * lim
            maxV = lim

    if cbarL != 0 and cbarM != 0:
        minV = np.int((np.nanmin(cbarL)))*1.0
        maxV = np.ceil(np.nanmax(cbarM))

    if ncolors == 0:
        cbar_levs = np.linspace(minV-ex, maxV+ex,(maxV-minV+2*ex)+1)
    else:
        cbar_levs = np.linspace(minV-ex,maxV+ex,ncolors*((maxV-minV+2*ex))+1.0)
        
    return cbar_levs

def shift_center(zfl_load,loncenter,latcenter,lon_arr,lat_arr):
    # first put lon,lat centered at 0 (originally runs 0 to 360)
    nlon,nlat=len(lon_arr),len(lat_arr)
    CENTER=np.zeros([nlon,nlat])*np.nan
    for i in range(0,len(lon_arr)):
        if i<nlon/2:
            CENTER[i+nlon/2,:]=zfl_load[i,:]
        if i>=nlon/2:
            CENTER[i-nlon/2,:]=zfl_load[i,:]
    plt_data=CENTER
    plt_lon=np.linspace(-180,180,nlon)
    
    if loncenter!=0:
        SHIFT=np.zeros([nlon,nlat])*np.nan
        plt_lon_n=np.zeros_like(plt_lon)

        nl=np.argmin(np.abs(plt_lon-loncenter))
        print plt_lon[nl]
        if nl<nlon/2:
            SHIFT[:(nlon/2-nl),:]=CENTER[(nlon-(nlon/2-nl)):,:]
            SHIFT[(nlon/2-nl):,:]=CENTER[:(nlon-(nlon/2-nl)),:]
            plt_lon_n[:(nlon/2-nl)]=plt_lon[(nlon-(nlon/2-nl)):]
            plt_lon_n[(nlon/2-nl):]=plt_lon[:(nlon-(nlon/2-nl))]
        if nl>nlon/2:
            SHIFT[:(3*nlon/2-nl),:]=CENTER[(nl-nlon/2):,:]
            SHIFT[(3*nlon/2-nl):,:]=CENTER[:(nl-nlon/2),:]
            plt_lon_n[:(3*nlon/2-nl)]=plt_lon[(nl-nlon/2):]
            plt_lon_n[(3*nlon/2-nl):]=plt_lon[:(nl-nlon/2)]

        plt_lon=plt_lon+loncenter
        plt_data=SHIFT
    return plt_data, plt_lon

def ortholine(R,lonl,latl,loncenter,latcenter):
    x=R*np.cos(latl)*np.sin(lonl-loncenter)
    y=R*(np.cos(latcenter)*np.sin(latl)-np.sin(latcenter)*np.cos(latl)*np.cos(lonl-loncenter))
    cosc=np.sin(latcenter)*np.sin(latl)+np.cos(latcenter)*np.cos(latl)*np.cos(lonl-loncenter)
    for i in range(0,len(cosc)):
        if cosc[i]<0:
            x[i]=np.nan
            y[i]=np.nan
    return x,y
    

def lat_v_lon_plot(x,y,z,cbar_levs,cmap_in,x_lab,y_lab,c_lab,savet,saven,bckgd):
    plt.figure(1,figsize=(10,6.25))
    plt.gcf().subplots_adjust(bottom=0.12,top=0.95,left=0.12,right=0.99)
    
    XX,YY=np.meshgrid(x,y)
        
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', cmap_in)
    p=plt.contourf(XX,YY,z.T,levels=cbar_levs,cmap=mymap,zorder=0)
    
    c=plt.colorbar(p)
    c.ax.tick_params(labelsize=18)

    plt.ylim(np.nanmin(y),np.nanmax(y))
    plt.xlim(np.nanmin(x),np.nanmax(x))

    plt.ylabel(y_lab,fontsize=20)
    plt.xlabel(x_lab,fontsize=20)
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)

    plt.figtext(0.84,0.5,c_lab,
                fontsize=20,rotation='vertical',ha='center',va='center')

    if savet==True:
        plt.savefig(saven,rasterized=True,transparent=bckgd)
        
    plt.show()
    
def lat_v_lon_ortho(x,y,z,loncenter,latcenter,cbar_levs,cmap_in,x_lab,y_lab,c_lab,savet,saven,bckgd):
    
    fig=plt.figure(1,figsize=(11,8))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])    
    ol=1 
    R=1.0
    
    lon_arr=x
    lat_arr=y
    
    nlon,nlat=len(lon_arr),len(lat_arr)
    
    #if units_a==1:  #convert to rads for ortho projection (needed to go to x/y)
    if np.nanmin(lat_arr)<-2.:
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
                plt_data0[i,j]=(z.T)[nlat+ol-i,j]
            if j>=nlon and i <nlat:
                plt_data0[i,j]=(z.T)[i,nlon+ol-j] 
            if j>=nlon and i>=nlat:
                plt_data0[i,j]=(z.T)[nlat+ol-i,nlon+ol-j] 
    plt_data0[:nlat,:nlon]=(z.T)

    MASK=np.ones_like(LON) # to mask points on unseen hemisphere

    X=R*np.cos(LAT)*np.sin(LON-loncenter)
    Y=R*(np.cos(latcenter)*np.sin(LAT)-np.sin(latcenter)*np.cos(LAT)*np.cos(LON-loncenter))
    cosc=np.sin(latcenter)*np.sin(LAT)+np.cos(latcenter)*np.cos(LAT)*np.cos(LON-loncenter)

    for i in range(0,X.shape[0]):  #checks for unseen hemisphere
        for j in range(0,X.shape[1]):
            if cosc[i,j]<0:
                MASK[i,j]=np.nan


    ax.axis('off')

    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', cmap_in)
    p=plt.contourf(X,Y,plt_data0*MASK,levels=cbar_levs,cmap=mymap,zorder=0,alpha=1.0)
    c=plt.colorbar(p)
    c.ax.tick_params(labelsize=18) 
    
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
    
    x,y=ortholine(R,120.*np.pi/180.,lat_arr,loncenter,latcenter)
    plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
    x,y=ortholine(R,-120.*np.pi/180.,lat_arr,loncenter,latcenter)
    plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
    
    x,y=ortholine(R,150.*np.pi/180.,lat_arr,loncenter,latcenter)
    plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
    x,y=ortholine(R,-150.*np.pi/180.,lat_arr,loncenter,latcenter)
    plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
    
    x,y=ortholine(R,180.*np.pi/180.,lat_arr,loncenter,latcenter)
    plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
    x,y=ortholine(R,-180.*np.pi/180.,lat_arr,loncenter,latcenter)
    plt.plot(x,y,linewidth=0.5,linestyle='--',color='white',alpha=0.7,zorder=10)
    
    
    plt.figtext(0.8,0.5,c_lab,fontsize=20,rotation='vertical',ha='center',va='center')
        
    
    if savet==True:
        plt.savefig(saven,rasterized=True,transparent=bckgd)
        
    plt.show()
    