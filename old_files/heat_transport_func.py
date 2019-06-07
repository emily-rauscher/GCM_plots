import numpy as np
from SurfacePressure_func import SurfPress
from altitude_func import altitude
import pickle

def meridtrans(temps,data50,sigma,latar,lonar,R,g,p0,oom):
    alt,temps=altitude(temps,data50,sigma,latar,lonar,R,g,p0,oom)
    cpt=0.286*R
    m=np.zeros_like(temps)*np.nan
    m=cpt*temps+g*alt
    return m,alt,temps
    
def viaeddies(temps,vwnd,surfp,sigma,latar,lonar,R,g,p0,oom):
    
    mtrn=meridtrans(temps,surfp,sigma,latar,lonar,R,g,p0,oom)[0]
    
    MEANm=np.nanmean(mtrn,axis=1)
    MEANv=np.nanmean(vwnd,axis=1)
   
    DEVm=np.zeros_like(mtrn)
    DEVv=np.zeros_like(vwnd)
    
    for j in range(0,mtrn.shape[1]):
        DEVm[:,j,:]=mtrn[:,j,:]-MEANm[:,:]
        DEVv[:,j,:]=vwnd[:,j,:]-MEANv[:,:]
        
    vedd=np.nanmean(DEVv*DEVm,axis=1)
    return vedd

def viameanfl(temps,vwnd,surfp,sigma,latar,lonar,R,g,p0,oom):
    mtrn=meridtrans(temps,surfp,sigma,latar,lonar,R,g,p0,oom)[0]
    
    MEANm=np.nanmean(mtrn,axis=1)
    MEANv=np.nanmean(vwnd,axis=1)
        
    vmnf=MEANv*MEANm
    return vmnf

def vertical_integral(data, sigma, p0, g, Rp, lat_arr):
    press=p0*sigma
    integrate=np.zeros_like(data)*np.nan
    lat_weight=2*np.pi*np.cos(lat_arr*np.pi/180.)*Rp**2.
    for i in range(0,len(sigma)):
        if i==0:
            integrate[i,:]=lat_weight*data[i,:]*(press[i])*(10**5.)/(g)*(1./Rp)
        else:
            integrate[i,:]=integrate[i-1,:]+lat_weight*data[i,:]*(press[i]-press[i-1])*(10**5.)/(g)*(1./Rp)
    return integrate
            
