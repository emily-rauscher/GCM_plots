from SurfacePressure_func import SurfPress
from load_data_noin import load_data
import numpy as np


def N2(path, runname, oom, surfp, grav,gascon, akap):
    data_50=SurfPress(path,runname,oom,surfp,False, False,'nan')
    runname,lon_arr,lat_arr,oom,surfp,p_BAR,data_26,data_lo,data_olr=load_data(path,runname,oom,surfp,False,False,False,'fort.26')

    #print data_50.shape
    #print data_26.shape
    
    Temp=data_26[:,:,:,5]
    
    #compute real pressure at each sigma level
    sigma=p_BAR/surfp
    Press=np.empty([len(lon_arr), len(lat_arr),len(p_BAR)])
    for l in range(0,len(p_BAR)):
        Press[:,:,l]=data_50[:,:,2]*sigma[l]
    
    #calc dtdp
    dtdp=np.empty(Temp.shape)*0.0
    N2=np.empty(Temp.shape)*0.0
    for i in range(0,len(lon_arr)):
        for j in range(0,len(lat_arr)):
            dtdp[:,i,j]=np.gradient(np.log(Temp[:,i,j]),np.log(Press[i,j,:]))
            N2[:,i,j]=grav**2./(gascon*(Temp[:,i,j]))*(akap-dtdp[:,i,j])
            #print N2[i,j,:]
            
    return lat_arr,lon_arr,N2
            