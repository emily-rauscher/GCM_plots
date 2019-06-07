import numpy as np
from SurfacePressure_func import SurfPress
import pickle

def altitude(temps,data50,sigma,latar,lonar,R,g,p0,oom):

    surfp=SurfPress(data50,oom,p0,False,False,'')[:,:,2]
    
    nlev,nlon,nlat=len(sigma),len(lonar),len(latar)
        
    TGR=np.nanmax(temps[:,:,nlev-1])
    z=np.empty([nlev,nlon,nlat])*np.nan
    z[nlev-1,:,:]=(R/g)*0.5*(temps[nlev-1,:,:]+TGR)*np.log(p0/surfp[:,:]/sigma[nlev-1])
    
    for i in range(nlev-2,-1,-1):
        z[i,:,:]=z[i+1,:,:]+(R/g)*0.5*(temps[i,:,:]+temps[i+1,:,:])*np.log(sigma[i+1]/sigma[i])
        
    return z,temps