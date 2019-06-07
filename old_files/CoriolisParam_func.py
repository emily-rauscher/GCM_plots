from load_data_noin import load_data
import numpy as np

def CoriolisParam(path,runname,oom,surfp,WW):
    runname,lon_arr,lat_arr,oom,surfp,p_BAR,data_26,data_lo,data_olr=load_data(path,runname,oom,surfp,False,False,False,'fort.26')
    
    lat_arr*=np.pi/180. #convert degrees to radians
    return lat_arr*180./np.pi,2.*WW*np.sin(lat_arr)