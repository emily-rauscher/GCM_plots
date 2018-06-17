import numpy as np
import math

from scipy.io import readsav

def load_data(path,LO,OLR,ver):
    filename=raw_input('Run Name?:')
    oom=input('Pressure OOM?:')
    surfp=input('Surface Press [bar]:')

    ############################################################

    data_26=readsav(path+filename+'/atm_wind_temp.sav')['xy']
    nlev,nlon,nlat,nparam=data_26.shape
    
    if LO==True:
        data_lo=readsav(path+filename+'/finalorb.sav')['all_dat']
    else:
        data_lo=np.empty(data_26.shape)
    if OLR==True:
        data_olr=readsav(path+filename+'/olr.sav')['olr']
    else:
        data_olr=np.empty(data_26.shape)

    if ver==True:
        print ' '
        print '--------------------------'
        print '|    ARRAY DIMENSIONS    |'
        print '--------------------------'
        print 'N_levels: ', nlev
        print 'N_lons:   ', nlon
        print 'N_lats:   ', nlat
        print 'N_params: ', nparam

    # nparam index: 
    #      0=lons
    #      1=lats
    #      2=levs
    #      3=u wind
    #      4=v wind
    #      5=temps

    sigma=np.empty([nlev])*0.0
    if oom==0:
        stp=1.0/(nlev+1.)
        sigma[nlev-1]=1.0-stp
        for n in range(nlev-2,-1,-1):
            sigma[n]=sigma[n+1]-stp
    if oom>0:
        stp=-1.0*oom/nlev
        sigma[nlev-1]=10.**(stp/2.)
        for n in range(nlev-2,-1,-1):
            sigma[n]=sigma[n+1]*10.**(stp)

    p_BAR=sigma*surfp
    if ver==True:
        print ' '
        print 'PRESSURE ARRAY: '
        print p_BAR  

    lat_arr=data_26[0,0,:,1]
    lon_arr=data_26[0,:,0,0]

    if ver==True:
        print ' '
        print 'LATITUDE ARRAY: '
        print lat_arr
        print ' '
        print 'LONGITUDE ARRAY: '
        print lon_arr  
    
    return filename,lon_arr,lat_arr,oom,surfp,p_BAR,data_26,data_lo,data_olr