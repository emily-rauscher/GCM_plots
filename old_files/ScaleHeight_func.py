from load_data_noin import load_data
import numpy as np

def ScaleHeight(path, runname, oom, surfp, grav,gascon):
    runname,lon_arr,lat_arr,oom,surfp,p_BAR,data_26,data_lo,data_olr=load_data(path,runname,oom,surfp,False,False,False,'fort.26')

    #print data_50.shape
    #print data_26.shape
    
    Temp=data_26[:,:,:,5]
    
    mu=(8.314/gascon)*1000.  #g/mol
    mh=1.6733*10**-24. #g
    grav=grav*100. #cm/s^2
    kb=1.38*10**-16  #cgs
    
    H=(kb*Temp/(mu*mh*grav))/(100.*1000.)
    return lat_arr,lon_arr,H