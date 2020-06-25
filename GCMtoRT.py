from astropy.io import ascii 
import numpy as np
import math
from scipy.io import readsav
from scipy.interpolate import griddata
from scipy import interpolate
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.font_manager import FontProperties
from astropy.table import Table,Column
from scipy.io import readsav

def altitude(path,runname,fort50,fort26,gasconst,grav,oom,tgr,nlay,surfp):

    #get surface pressures
    with open(path+runname+'/'+fort50,'r') as data_50:  #fort_50 long1 lat1 pressure 
        specificp=np.zeros((48,96))  #lat,long               #long1 lat2 pressure
        acount=0
        bcount=0
        next(data_50)
        for line in data_50:
            p=line.split()
            if bcount<48 and acount<96:
                specificp[bcount][acount]=((float(p[2]))+1)*surfp
                bcount=bcount+1
            else:
                if acount<96:
                    acount=acount+1
                    bcount=0
                    if acount<96:
                        specificp[bcount][acount]=((float(p[2]))+1)*surfp
                        bcount=bcount+1
                
    #get T
    with open(path+runname+'/'+fort26) as f:
        first_line=f.readline()
        nlat,nlon,nlev=first_line.split()
        nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)

    f.close()
    data26=np.empty([nlon*nlat*nlev, 6])
    
    l=0
    lp=0
    with open(path+runname+'/'+fort26) as f:
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
                print ('       END OF FILE: DONE')
                break
            l+=1
    f.close()

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

    data_26=np.empty([nlev,nlon,nlat,6])
    for l in range(0,data26.shape[0]):
        lon,lat,lev=data26[l,:3]
        lon_i,lat_i,lev_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0],np.where(lev_arr==lev)[0][0]
        data_26[lev_i,lon_i,lat_i,:]=data26[l,:]
    ############################################################
    nlev,nlon,nlat,nparam=data_26.shape
    temps=data_26[:,:,:,5]
    
    
    #set sigma
    sigma=np.empty([nlay])*0.0
    if oom>0: #setting up pressure values 
        stp=-1.0*oom/nlay
        sigma[nlay-1]=10.**(stp/2.)
        for n in range(nlay-2,-1,-1):
            sigma[n]=sigma[n+1]*10.**(stp)
    p_BAR=sigma*surfp
    
    #create array to hold heights
    z=np.zeros((nlon,nlat,nlay))
    #set altitude of the first level (up from base=p0, where T=TGR)
    z[:,:,nlay-1]=(gasconst/grav) *.5*(temps[nlay-1,:,:] + tgr) * np.log(surfp/specificp.T/sigma[nlay-1])
    #integrate hydrostatic to solve for higher levels
    start=nlev-2
    while start >= 0:
        z[:,:,start]=z[:,:,start+1] + (gasconst/grav) *0.5*(temps[start,:,:]+temps[start+1,:,:]) * np.log(sigma[start+1]/sigma[start])
        start=start - 1
      #making vertical.txt  
    vertinfo=np.zeros((nlay*nlat*nlon,3))
    h= open(path+runname+'/vertical.txt', 'w')   
    h.write("%5.4f %5.4f %5.4f \n" % (nlat,nlon,nlay))
    for ii in range(nlay):
       # print(ii)
        for jj in range(nlon):
            for kk in range(nlat):
                h.write("%10.5e %10.5e %10.5e \n" % (sigma[ii], sigma[ii]*specificp[kk,jj], z[jj,kk,ii]))
    h.write("\n")
    h.write(" R,g,TGR,p0,OOM: %5.2f %5.4f %5.2f %5.2f %5.2f" % (gasconst, grav, tgr, surfp, oom))
    h.close()
    return data_26,vertinfo

def readtp(path,runname):
    for phasenumb in range(90):
        if phasenumb<10: #file 5 is reported at 2605
            ident='0'+ str(phasenumb)
        else:
            ident=str(phasenumb)
        for26=str(path+runname+'/fort.26'+ident)
        vertic=str(path+runname+'/vertical.txt')
        
        #getting the info from fort.26
        with open(for26) as f:
            first_line=f.readline()
            nlat,nlon,nlev=first_line.split()
            nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)

            f.close()
        data26=np.empty([nlon*nlat*nlev, 6])
    
        l=0
        lp=0
        with open(for26) as f:
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
                    #print ('       END OF FILE: DONE')
                    break
                l+=1
        f.close()

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

        data_26=np.empty([nlev,nlon,nlat,6])
        for l in range(0,data26.shape[0]):
            lon,lat,lev=data26[l,:3]
            lon_i,lat_i,lev_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0],np.where(lev_arr==lev)[0][0]
            data_26[lev_i,lon_i,lat_i,:]=data26[l,:]
            
        #getting vertical data
        ver=pd.read_csv(vertic, delim_whitespace=True, header=0,skipfooter=1,engine='python')
        #print (np.shape(ver)) #196460,3 
        sigma=ver.values[:,0]
        p=ver.values[:,1]
        alt=ver.values[:,2]
        nlev,nlon,nlat,nparam=data_26.shape
        temps=data_26[:,:,:,5]
        #print (np.shape(temps))
        #get rid of 0 entries in t_temp Confused about this though,
        for aa in range(nlev):
            for bb in range(nlon):
                for cc in range(nlat):
                    if temps[aa,bb,cc]== 0:
                        print ('got a 0')
                        temps[aa,bb,cc]==1
                        
        ntau=60
        level_ll=np.linspace(0,44,45) + 1
        lat_ll=np.zeros((nlat))
        long_ll=np.zeros((nlon+1))
        alt_ll=np.zeros((nlat, nlon, nlev))
        p_ll=np.zeros((nlat,nlon,nlev))
        t_ll=np.zeros((nlat,nlon,nlev))
        ew_vel_ll=np.zeros((nlat, nlon, nlev))
        ns_vel_ll=np.zeros((nlat, nlon, nlev))
        n=np.copy(nlat)
        
        for dd in range(nlat):
            lat_ll[dd]=lat_arr[nlat-1-dd]
        lon_ll=np.copy(lon_arr)
        
        break   
        

def createsobslon(path,runname,prot,porb): #goes through all fort.26xx and returns 3 arrays
    phasenumb=0                            #1: sub-observer longitudes, sub-stellar longitudes, phase values
    sobslon=np.zeros((90))
    dayarr=np.zeros((90))
    sslon=np.zeros((90))
    dayphase=np.zeros((90))
    while phasenumb<90:
        if phasenumb<10: #file 5 is reported at 2605
            ident='0'+ str(phasenumb)
        else:
            ident=str(phasenumb)
        with open(path+runname+'/fort.26'+ident) as g: 
            for line in g:
                pass
            last = line
            sslon[phasenumb]= (''.join(last[49:57])) #grabbing the sub-stellar longitudes
            dayarr[phasenumb]=(''.join(last[18:27])) #grabbing the day outputs
        g.close()
        dayphase[phasenumb]=(float(dayarr[phasenumb]-dayarr[0]))*(prot/porb) #converting days to phase with 0 being the first file       
        sobslon[phasenumb]=sslon[phasenumb]+180.-(360.*dayphase[phasenumb])  # sub-observer point moves W as planet spins E
        sobslon[phasenumb] = (((sobslon[phasenumb] % 360)+360.) % 360)
        phasenumb=phasenumb+1
        #if phasenumb % 10 ==0:
           # print (phasenumb)
    return sobslon,sslon,dayphase


def avginterpolate(path,runname,prot,porb,phasenumb,rotvalue,nlon,tdsync,name):
    if phasenumb<10: #file 5 is reported at 2605
        ident='0'+ str(phasenumb)
    else:
        ident=str(phasenumb)
     
    file_orig = np.loadtxt(path+runname+'/tprofiles/t_p_3D_'+str(phasenumb)+'.dat') #imports original data file 
    file_orig = np.asarray(file_orig)
    file_orig = file_orig.astype(np.float)
    if tdsync== True:  
        file_sync=np.copy(file_orig)
    else:
        file_sync = np.loadtxt(path+'test03/tprofiles/t_p_3D_'+str(phasenumb)+'.dat')  #need to read in sync case in order to get lat and lon grid 
    file_sync = np.asarray(file_sync)
    file_sync = file_sync.astype(np.float)
    nlat = 48  
    level_length = nlat*nlon 
    columns = 4 #number of columns of physical parameters beginning after altitude column 
    levels = 60 #number of altitude levels
    alt_col = 4 #column number for altitude -- parameter columns begin after this
    porbrat=porb/prot  #correct ratio (check fort.7 if confused)

    latarrayrad = file_orig[:,0]*(np.pi/180)
    lonarrayrad = file_orig[:,1]*(np.pi/180)
    filelength = len(lonarrayrad) 
    dayphase=np.zeros((90))
    lat_prerot = np.zeros([nlat,nlon])
    lon_prerot = np.zeros([nlat,nlon])
    rot_lat = np.zeros([level_length,])
    rot_lon = np.zeros([level_length,])
    orig_lat = np.zeros([level_length,])
    orig_lon = np.zeros([level_length,])
    param_out = np.zeros([48*97, columns]) #used to be level_length
    

    for ilev in range (levels):
        f = open(path+runname+'/interfile/nophase_level{}.txt'.format(1+ilev), 'w')
        level_check = file_orig[:,2][0+ilev:filelength:levels]
        altitude_sync = file_sync[:,3][0+ilev:48*97*60:levels]
        level_sync = file_sync[:,2][0+ilev:48*97*60:levels]
        altitude = file_orig[:,3][0+ilev:filelength:levels] 
        rot_lat = file_orig[:,0][0+ilev:filelength:levels]
        rot_lon = file_orig[:,1][0+ilev:filelength:levels]
        orig_lat = file_orig[:,0][0+ilev:filelength:levels]
        sync_lat = file_sync[:,0][0+ilev:48*97*60:levels]
        sync_lon = file_sync[:,1][0+ilev:48*97*60:levels]
        file_lon = file_orig[:,1][0+ilev:filelength:levels] 
        orig_lon = file_orig[:,1][0+ilev:filelength:levels] #print original lon range back to file
        
        ##Rotating the longitude grid####
        rot_lon=orig_lon+(180-rotvalue) #longitude rotation by sub observer longitude 
        check=0
        while check< len(rot_lon):   #need to maybe make sure 0 and 360 are covered in the regular lons
            if rot_lon[check]<0:
                rot_lon[check]=360+rot_lon[check]
            if rot_lon[check]>360:
                rot_lon[check]=rot_lon[check]-360     
            check=check+1
        #copy minimum  longitude but rename it 360 +minimum for interpolation purposes 
        small_lon=min(rot_lon)
        big_lon=max(rot_lon)
        lon_count=0
        z=0
        lon_addition=np.zeros((2*nlat))
        lat_addition=np.zeros((2*nlat))
        endpoint_orig=np.zeros((2*nlat,8))
        rotlon_extended=[]
        rotlat_extended=[]
        while lon_count< len(rot_lon):
            if rot_lon[lon_count]==small_lon:
                zero_row=np.copy(file_orig[ilev+lon_count*levels,:]) #need to find index of this copy it from file orig and then append to file orig
                zero_row[1]=360+small_lon #need to add it to file_orig, rot_lon and rot_lat
                lon_addition[z]=(zero_row[1])
                lat_addition[z]=(zero_row[0])
                endpoint_orig[z,:]=(zero_row)
                z=z+1
            if rot_lon[lon_count]==big_lon:
                neg_row=np.copy(file_orig[ilev+lon_count*levels,:])
                neg_row[1]=big_lon-360
                lon_addition[z]=(neg_row[1])
                lat_addition[z]=(neg_row[0])
                endpoint_orig[z,:]=(neg_row)
                z=z+1
            lon_count=lon_count+1
        rotlon_extended=np.concatenate((rot_lon,lon_addition))
        rotlat_extended=np.concatenate((rot_lat,lat_addition))  
        for icol in range(columns):
            parameter_orig = file_orig[:,alt_col+icol][0+ilev:filelength:levels] 
            params=[]
            parameter_full=[]
            params=endpoint_orig[:,alt_col+icol]
            parameter_full=np.concatenate((parameter_orig,params))
            int_loop = griddata((rotlon_extended,rotlat_extended),parameter_full,(sync_lon,sync_lat),method='cubic',fill_value=parameter_orig[30]) #the correct interp
            param_out[:,icol] = int_loop
    #writing to file 
        for i in range(48*97): #used to be level length but change for sync case 
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(format(sync_lat[i], 'E'), format(sync_lon[i], 'E'), format(level_sync[i]), format(altitude_sync[i], 'E'), format(param_out[i,0], 'E'), format(param_out[i,1], 'E'), format(param_out[i,2], 'E'), format(param_out[i,3], 'E')))
        f.close()
    final = open(path+runname+'/interfile/'+runname+ident+'inter.dat', 'w')
    for r in range(48*97):
        for j in range (levels):
            lev = open(path+runname+'/interfile/nophase_level{}.txt'.format(1+j), 'r')
            levlines = list(lev)
            final.write(levlines[r])
    final.close()
    print ('interpolation finished')
    old_file=path+runname+'/interfile/'+runname+ident+'inter.dat'
    data = np.loadtxt(old_file) 
    x=np.shape(data)[0]
    data = ascii.read(old_file, names=('lat', 'lon', 'level', 'alt', 'pressure', 'temp', 'ew', 'ns'))
    #adding zeros for  13 cloud columns
    clouds1=Column(np.zeros(x),name='cl1')
    clouds2=Column(np.zeros(x),name='cl2')
    clouds3=Column(np.zeros(x),name='cl3')
    clouds4=Column(np.zeros(x),name='cl4')
    clouds5=Column(np.zeros(x),name='cl5')
    clouds6=Column(np.zeros(x),name='cl6')
    clouds7=Column(np.zeros(x),name='cl7')
    clouds8=Column(np.zeros(x),name='cl8')
    clouds9=Column(np.zeros(x),name='cl9')
    clouds10=Column(np.zeros(x),name='cl10')
    clouds11=Column(np.zeros(x),name='cl11')
    clouds12=Column(np.zeros(x),name='cl12')
    clouds13=Column(np.zeros(x),name='cl13')
    data.add_column(clouds1)
    data.add_column(clouds2)
    data.add_column(clouds3)
    data.add_column(clouds4)
    data.add_column(clouds5)
    data.add_column(clouds6)
    data.add_column(clouds7)
    data.add_column(clouds8)
    data.add_column(clouds9)
    data.add_column(clouds10)
    data.add_column(clouds11)
    data.add_column(clouds12)
    data.add_column(clouds13)

  
    #remove rows where longitude is 360, if these exist in file 

    keep = np.where(data['lon'] != 360.0) #keep data where lon is not 360

    good_data = data[keep]

    #make array from table in same shape as original data 

    data_array = np.array([good_data['lat'],good_data['lon'],good_data['level'],good_data['alt'],good_data['pressure'], good_data['temp'],good_data['ew'],good_data['ns'],good_data['cl1'], good_data['cl2'], good_data['cl3'], good_data['cl4'], good_data['cl5'], good_data['cl6'], good_data['cl7'], good_data['cl8'], good_data['cl9'], good_data['cl10'], good_data['cl11'], good_data['cl2'], good_data['cl3'] ])
    best_data=(np.transpose(data_array))


    #copy data and add 720 values 
    data2 = np.copy(best_data)
    data3 = []

    for row in data2:
        row[1] += 360.
        if row[1] == 360.:
            last_row = np.copy(row)
            last_row[1] += 360.
            data3.append(last_row)
        

    double_data = np.concatenate((best_data,data2,data3)) #add data3 back?

    #sort data in order of ascending latitude

    double_data = double_data[double_data[:,0].argsort()]
    # sort data to match original file format (this might not be the most efficient way to do this)
    nlat = 48
    nlon = 195 #used to be nlon*2 +1 
    for i in range(nlat):
        chunk = double_data[i*len(double_data)/nlat:(i+1)*len(double_data)/nlat]
        chunk = chunk[chunk[:,1].argsort()]
        for j in range(nlon):
            sub_chunk = chunk[j*len(chunk)/nlon:(j+1)*len(chunk)/nlon]
            sub_chunk = sub_chunk[sub_chunk[:,2].argsort()]
            chunk[j*len(chunk)/(nlon):(j+1)*len(chunk)/(nlon)] = sub_chunk    
        double_data[i*len(double_data)/nlat:(i+1)*len(double_data)/nlat] = chunk
    
    # save the new file -- format will need to be adjusted, depending on original file
    np.savetxt(path+runname+'/t_p_3D_'+runname+'rotated'+str(name)+'.dat', double_data,           fmt='%.5E\t%.5E\t%i\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t%.5E\t')


    return final