import numpy as np
import math

from scipy.io import readsav

def load_data(path,LO,OLR,ver):
    runname=raw_input('Run Name?:')
    oom=input('Pressure OOM?:')
    surfp=input('Surface Press [bar]:')

    ####################### READ 26 #########################
    with open(path+runname+'/fort.26') as f:
        first_line=f.readline()
        nlat,nlon,nlev=first_line.split()
        nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)
        print '  '
        print ' ....reading fort.26'
        print '       nlat=', nlat, 'nlon=', nlon, 'nlev=', nlev
    f.close()
    
    data26=np.empty([nlon*nlat*nlev, 6])
    
    l=0
    lp=0
    with open(path+runname+'/fort.26') as f:
        for line in f:
            if l==0:
                l+=1
                continue
            elif l%2==1 and l<=nlon*nlat*nlev*2.:
                line_pair=np.empty([6])
                lon, lat, lev, u, v = line.split()
                line_pair[:5] = np.float32(lon), np.float32(lat), int(lev), np.float32(u), np.float32(v)
            elif l%2==0 and l<=nlon*nlat*nlev*2.:
                line_pair[5]=np.float32(line)
                data26[lp,:]=line_pair
                lp+=1
            elif l>nlon*nlat*nlev*2.:
                print '       END OF FILE: DONE'
                break
            l+=1

    lon_arr_f=data26[:,0]
    lon_arr=np.array([])
    for l in range(0,len(lon_arr_f)):
        el=lon_arr_f[l]
        if not el in lon_arr:
            lon_arr=np.append(lon_arr,el)

    lat_arr_f=data26[:,1]
    lat_arr=np.array([])
    for l in range(0,len(lat_arr_f)):
        el=lat_arr_f[l]
        if not el in lat_arr:
            lat_arr=np.append(lat_arr,el)

    lev_arr_f=data26[:,2]
    lev_arr=np.array([])
    for l in range(0,len(lev_arr_f)):
        el=lev_arr_f[l]
        if not el in lev_arr:
            lev_arr=np.append(lev_arr,el)

    #data26=np.reshape(data26,[nlon,nlat,nlev,6])
    #print data26.shape
    data_26=np.empty([nlev,nlon,nlat,6])
    for l in range(0,data26.shape[0]):
        lon,lat,lev=data26[l,:3]
        lon_i,lat_i,lev_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0],np.where(lev_arr==lev)[0][0]
        data_26[lev_i,lon_i,lat_i,:]=data26[l,:]
    ############################################################
    #data_26=readsav(path+runname+'/atm_wind_temp.sav')['xy']
    nlev,nlon,nlat,nparam=data_26.shape
    
    if LO==True:
        data_lo=readsav(path+runname+'/finalorb.sav')['all_dat']
    else:
        data_lo=np.empty(data_26.shape)
    if OLR==True:
        data_olr=readsav(path+runname+'/olr.sav')['olr']
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
    
    return runname,lon_arr,lat_arr,oom,surfp,p_BAR,data_26,data_lo,data_olr