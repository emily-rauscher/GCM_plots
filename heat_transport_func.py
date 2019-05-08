import numpy as np
from SurfacePressure_func import SurfPress
from altitude_func import altitude
import pickle

def meridtrans(path,runpath,lo,R,g,p0,oom):
    alt,temps=altitude(path,runpath,lo,R,g,p0,oom)
    cpt=0.286*R
    m=np.zeros_like(temps)*np.nan
    m=cpt*temps+g*alt
    return m,alt,temps
    
def viaeddies(path,runpath,lo,R,g,p0,oom):
    tp=path+runpath
    if lo==True:
        vwnd=pickle.load(open(tp+'/VWD_lo.txt', 'rb'))
    if lo==False:
        vwnd=pickle.load(open(tp+'/fort26.txt', 'rb'))[:,:,:,4]
    mtrn=meridtrans(path,runpath,lo,R,g,p0,oom)[0]
    
    MEANm=np.nanmean(mtrn,axis=1)
    MEANv=np.nanmean(vwnd,axis=1)
   
    DEVm=np.zeros_like(mtrn)
    DEVv=np.zeros_like(vwnd)
    
    for j in range(0,mtrn.shape[1]):
        DEVm[:,j,:]=mtrn[:,j,:]-MEANm[:,:]
        DEVv[:,j,:]=vwnd[:,j,:]-MEANv[:,:]
        
    vedd=np.nanmean(DEVv*DEVm,axis=1)
    return vedd

def viameanfl(path,runpath,lo,R,g,p0,oom):
    tp=path+runpath
    if lo==True:
        vwnd=pickle.load(open(tp+'/VWD_lo.txt', 'rb'))
    if lo==False:
        vwnd=pickle.load(open(tp+'/fort26.txt', 'rb'))[:,:,:,4]
    mtrn=meridtrans(path,runpath,lo,R,g,p0,oom)[0]
    
    MEANm=np.nanmean(mtrn,axis=1)
    MEANv=np.nanmean(vwnd,axis=1)
        
    vmnf=MEANv*MEANm
    return vmnf
