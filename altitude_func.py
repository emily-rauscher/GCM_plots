import numpy as np
from SurfacePressure_func import SurfPress
import pickle

def altitude(path,runpath,lo,R,g,p0,oom):
    tp=path+runpath
    sigma=pickle.load(open(tp+'/pres_lon_lat.txt', 'rb'))[0]/p0
    latar=pickle.load(open(tp+'/pres_lon_lat.txt', 'rb'))[2]
    lonar=pickle.load(open(tp+'/pres_lon_lat.txt', 'rb'))[1]
    nlat,nlon,nlev=len(latar),len(lonar),len(sigma)
    if lo==True:
        temps=pickle.load(open(tp+'/TMP_lo.txt', 'rb'))
        surfp=SurfPress(path,runpath,oom,p0,False,False,'')[:,:,2] #update last orb to do fort.50 too
    if lo==False:
        temps=pickle.load(open(tp+'/fort26.txt', 'rb'))[:,:,:,5]
        surfp=SurfPress(path,runpath,oom,p0,False,False,'')[:,:,2]
        
    TGR=np.nanmax(temps[:,:,nlev-1])
    z=np.empty([nlev,nlon,nlat])*np.nan
    z[nlev-1,:,:]=(R/g)*0.5*(temps[nlev-1,:,:]+TGR)*np.log(p0/surfp[:,:]/sigma[nlev-1])
    
    for i in range(nlev-2,-1,-1):
        z[i,:,:]=z[i+1,:,:]+(R/g)*0.5*(temps[i,:,:]+temps[i+1,:,:])*np.log(sigma[i+1]/sigma[i])
        
    return z,temps