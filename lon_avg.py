import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from matplotlib.font_manager import FontProperties
font0=FontProperties()
font=font0.copy()
font.set_family('serif')

fontb=font.copy()
fontb.set_weight('bold')

params = {'font.family': 'serif',}
matplotlib.rcParams.update(params)

from load_data import load_data

#data_26 index: Nlev,nlon,nlat,nparam,

# nparam index: 
#      0=lons
#      1=lats
#      2=levs
#      3=u wind
#      4=v wind
#      5=temps
#def lon_avg(plot,lev,lo,data_t21,data_t42,data_t63,lats_t21,lats_t42,lats_t63):
def lon_avg_comp(path,plot,lev,lo,ln,tn):
    nnames=0
    print '********** Maximum of 4 RUNS currently **********'
    print '*************************************************'
    ###################
    enter='y'#raw_input('Enter another run name?? (y or n)')
    if 'y' in enter:
        runname_1,lon_arr_1,lat_arr_1,oom_1,surfp_1,p_BAR_1,data_26_1,data_lo_1=load_data(path,lo,False)
        nnames+=1
        ###################
        print '*************************************************'
        enter=raw_input('Enter another run name?? (y or n)')
        if 'y' in enter:
            runname_2,lon_arr_2,lat_arr_2,oom_2,surfp_2,p_BAR_2,data_26_2,data_lo_2=load_data(path,lo,False)
            nnames+=1
            print '*************************************************'
            ###################
            enter=raw_input('Enter another run name?? (y or n)')
            if 'y' in enter:
                runname_3,lon_arr_3,lat_arr_3,oom_3,surfp_3,p_BAR_3,data_26_3,data_lo_3=load_data(path,lo,False)
                nnames+=1
                print '*************************************************'
                ###################
                enter=raw_input('Enter another run name?? (y or n)')
                if 'y' in enter:
                    runname_4,lon_arr_4,lat_arr_4,oom_4,surfp_4,p_BAR_4,data_26_4,data_lo_4=load_data(path,lo,False)
                    nnames+=1
                    print '*************************************************'
                else:
                    print 'DONE WITH RUNS. ENTERED ',nnames,' RUN NAMES'
            else:
                print 'DONE WITH RUNS. ENTERED ',nnames,' RUN NAMES'
        else:
            print 'DONE WITH RUNS. ENTERED ',nnames,' RUN NAMES'
    else:
        print 'DONE WITH RUNS. ENTERED ',nnames,' RUN NAMES'

    #########################    
    if plot==0: #TEMPERATURE
        if lo==True:
            ind=0
        else:  
            ind=5
    if plot==1:
        if lo==True:
            ind=1
        else:
            ind=3 # uwind
    if plot==2:
        if lo==True:
            ind=2
        else:
            ind=4 # vwind
    
    if lo==True:
        data_1=np.median((np.nanmedian(data_lo_1,axis=4))[lev,:,:,ind],axis=0)
        if nnames>1:
            data_2=np.median((np.nanmedian(data_lo_2,axis=4))[lev,:,:,ind],axis=0)
            if nnames>2:
                data_3=np.median((np.nanmedian(data_lo_3,axis=4))[lev,:,:,ind],axis=0)
                if nnames>3:
                    data_4=np.median((np.nanmedian(data_lo_4,axis=4))[lev,:,:,ind],axis=0)
    
    else:
        data_1=np.median(np.copy(data_26_1[lev,:,:,ind]),axis=0)
        if nnames>1:
            data_2=np.median(np.copy(data_26_2[lev,:,:,ind]),axis=0)
            if nnames>2:
                data_3=np.median(np.copy(data_26_3[lev,:,:,ind]),axis=0)
                if nnames>3:
                    data_4=np.median(np.copy(data_26_4[lev,:,:,ind]),axis=0)
    
    
    plt.figure(figsize=(6.25,10))
    plt.gcf().subplots_adjust(bottom=0.08,top=0.97,left=0.17,right=0.97)
    
    ncolor=40
    color_list = plt.cm.plasma(np.linspace(0., 1, ncolor))
    if nnames==1:
        color1=color_list[int(ncolor/2)]
    if nnames==2:
        color1=color_list[int(ncolor/3)]
        color2=color_list[int(2.*ncolor/3)]
    if nnames==3:
        color1=color_list[int(8*ncolor/10)]
        color2=color_list[int(ncolor/2)]
        color3=color_list[int(ncolor/10)]
    if nnames==4:
        color1=color_list[int(8*ncolor/10)]
        color2=color_list[int(6*ncolor/10)]
        color3=color_list[int(3*ncolor/10)]
        color4=color_list[int(1*ncolor/10)]
    
    tn=tn
    ln=ln
    plt.plot(data_1,lat_arr_1,linewidth=5.0,linestyle=':',color=color1)
    plt.figtext(ln,tn,runname_1,fontsize=20,fontproperties=fontb,color=color1)
    if nnames>1:
        plt.plot(data_2,lat_arr_2,linewidth=5.0,linestyle='--',color=color2)
        plt.figtext(ln,tn-0.04,runname_2,fontsize=20,fontproperties=fontb,color=color2)
        if nnames>2:
            plt.plot(data_3,lat_arr_3,linewidth=5.0,linestyle='-',color=color3)
            plt.figtext(ln,tn-0.04*2.,runname_3,fontsize=20,fontproperties=fontb,color=color3)
            if nnames>3:
                plt.plot(data_4,lat_arr_4,linewidth=5.0,linestyle='-.',color=color4)
                plt.figtext(ln,tn-0.04*3.,runname_4,fontsize=20,fontproperties=fontb,color=color4)
    
    plt.ylabel('Latitude [degrees]',fontsize=20)
    if plot==0:
        plt.xlabel('Temp [K]',fontsize=20)
    if plot==1:
        plt.xlabel('E-W Wind [m/s]',fontsize=20)
    if plot==2:
        plt.xlabel('N-S Wind [m/s]',fontsize=20)
    
    plt.yticks(fontsize=18,fontproperties=font)
    plt.xticks(fontsize=18,fontproperties=font)
    
    if plot==0:
        plt.savefig(path+runname_1+'/LongAvg_Temps_Comparison.pdf')
        if nnames>1:
            plt.savefig(path+runname_2+'/LongAvg_Temps_Comparison.pdf')
            if nnames>2:
                plt.savefig(path+runname_3+'/LongAvg_Temps_Comparison.pdf')
                if nnames>3:
                    plt.savefig(path+runname_4+'/LongAvg_Temps_Comparison.pdf')
    if plot==1:
        plt.savefig(path+runname_1+'/LongAvg_UWinds_Comparison.pdf')
        if nnames>1:
            plt.savefig(path+runname_2+'/LongAvg_UWinds_Comparison.pdf')
            if nnames>2:
                plt.savefig(path+runname_3+'/LongAvg_UWinds_Comparison.pdf')
                if nnames>3:
                    plt.savefig(path+runname_4+'/LongAvg_UWinds_Comparison.pdf')
    if plot==2:
        plt.savefig(path+runname_1+'/LongAvg_VWinds_Comparison.pdf')
        if nnames>1:
            plt.savefig(path+runname_2+'/LongAvg_VWinds_Comparison.pdf')
            if nnames>2:
                plt.savefig(path+runname_3+'/LongAvg_VWinds_Comparison.pdf')
                if nnames>3:
                    plt.savefig(path+runname_4+'/LongAvg_VWinds_Comparison.pdf')
    
    plt.show()

def lon_avg(plot,path,data,lon_arr,lat_arr,lev,lo,ln,tn):
    if lo==True:
        data_lo_1=data
    else:
        data_26_1=data
    lon_arr_1=lon_arr
    lat_arr_1=lat_arr
    #########################    
    if plot==0: #TEMPERATURE
        if lo==True:
            ind=0
        else:  
            ind=5
    elif plot==1:
        if lo==True:
            ind=1
        else:
            ind=3 # uwind
    elif plot==2:
        if lo==True:
            ind=2
        else:
            ind=4 # vwind
    
    if lo==True:
        data_1=np.median((np.nanmedian(data_lo_1,axis=4))[lev,:,:,ind],axis=0)
    else:
        print lev,ind
        data_1=np.median(np.copy(data_26_1[lev,:,:,ind]),axis=0)
    
    
    plt.figure(figsize=(6.25,8))
    plt.gcf().subplots_adjust(bottom=0.08,top=0.97,left=0.17,right=0.97)
    
    ncolor=40
    color_list = plt.cm.plasma(np.linspace(0., 1, ncolor))
    color1=color_list[int(ncolor/2)]

    
    tn=tn
    ln=ln
    plt.plot(data_1,lat_arr_1,linewidth=5.0,linestyle='-',color=color1)

    
    plt.ylabel('Latitude [degrees]',fontsize=20)
    if plot==0:
        plt.xlabel('Temperature [K]',fontsize=20)
    if plot==1:
        plt.xlabel('E-W Wind [m/s]',fontsize=20)
    if plot==2:
        plt.xlabel('N-S Wind [m/s]',fontsize=20)
    
    plt.yticks(fontsize=18,fontproperties=font)
    plt.xticks(fontsize=18,fontproperties=font)
    
    if plot==0:
        plt.savefig(path+'/LongAvg_Temps.pdf')
    if plot==1:
        plt.savefig(path+'/LongAvg_UWinds.pdf')
    if plot==2:
        plt.savefig(path+'/LongAvg_VWinds.pdf')

    plt.show()

