import numpy as np
import math
from scipy.interpolate import griddata

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches

from matplotlib.font_manager import FontProperties
font0=FontProperties()
font=font0.copy()
font.set_family('serif')

fontb=font.copy()
fontb.set_weight('bold')

params = {'font.family': 'serif',}
matplotlib.rcParams.update(params)


def KEplot(path,runname,oom,press,ndays,lay,savename,minv,maxv):
    filename=runname
    oom=oom#input('Pressure OOM?:')
    surfp=press#input('Surface Press [bar]:')
    nlay=lay
    nlev=lay
    minV=minv
    maxV=maxv
    nday=ndays#input('Number of days:')
    with open(path+filename+'/'+'fort.52','r') as data_52:
        ke=[] #kinetic energy 
        cpt=np.zeros((nlay,nday+91))*np.nan #kinetic energy and cpt
        dayval=[] #column corresponding to the day 
        daycount=0
        laycount=0
        for line in data_52:
            p=line.split()
            #print laycount,daycount
            if laycount < nlay-1: 
                ke.append(float(p[0]))
                if daycount != nday+91:
                    cpt[laycount,daycount]=(float(p[0]))
                    dayval.append(daycount)
                    laycount=laycount+1
            else:
                laycount=0
                ke.append(float(p[0]))
                if daycount != nday+91:
                    cpt[laycount,daycount]=(float(p[0]))
                    dayval.append(daycount)
                    daycount=daycount+1

        cpt=cpt*surfp*100000
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
    cpt=cpt+.01
    cpt=np.log(cpt)
    daylist=np.arange(3,nday+1,1)
    #print len(daylist)
    
    newke=cpt[:,3:nday+1]
    #print newke
    
    #messing with the color bar
    colors1 = plt.cm.YlGnBu_r(np.linspace(0, 1, 128))
    colors2 = plt.cm.YlOrBr(np.linspace(0., 1, 128))
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
   
    ex=1
    ncolors=0
   
    
    
    plt.figure(figsize=(8.25,8))
    daylist=daylist[:]
    p_BAR=p_BAR[:nlay]
    #p_BAR=p_BAR[::-1]
    day,PRESS_P=np.meshgrid(daylist,p_BAR)
    
    #putting lower levels at bottom of graph
    #newke=np.flip(newke,0)
    
    if minV!=0 and maxV!=0:
        cbar_levs=np.linspace(minV-ex,maxV+ex,(maxV-minV+2*ex)+1)
        p=plt.contourf(day,PRESS_P,newke,levels=cbar_levs,cmap=mymap)
    else:
        cbar_levs=np.linspace(np.round(np.nanmin(newke),2)-1,np.round(np.nanmax(newke),2)+1,100)
    p=plt.contourf(day,PRESS_P,newke,levels=cbar_levs,cmap=mymap)
    c=plt.colorbar(p)
    c.ax.tick_params(labelsize=18)
    plt.grid(color='white',linewidth=0.5,linestyle='--',alpha=0.5, zorder=10)
    if oom>0:
        plt.yscale('log')
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure [bar]', fontsize=20)
    plt.xlabel('Planet Day', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    plt.savefig(savename,rasterized=True)
    plt.show()
    plt.close()
                       
            