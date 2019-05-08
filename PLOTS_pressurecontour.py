import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
#%matplotlib inline
from matplotlib.font_manager import FontProperties

import os
from matplotlib import font_manager as fm, rcParams


apath_r=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/AmaticSC-Regular.ttf')
apath_b=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/Amatic-Bold.ttf')
aprop_r=fm.FontProperties(fname=apath_r)
aprop_b=fm.FontProperties(fname=apath_b)

jpath_r=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/JosefinSlab-Regular.ttf')
jpath_sb=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/JosefinSlab-SemiBold.ttf')
jpath_b=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/JosefinSlab-Bold.ttf')
jprop_r=fm.FontProperties(fname=jpath_r)
jprop_sb=fm.FontProperties(fname=jpath_sb)
jprop_b=fm.FontProperties(fname=jpath_b)

opath_r=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/Oswald-Regular.ttf')
opath_sb=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/Oswald-DemiBold.ttf')
opath_b=os.path.join(rcParams["datapath"],'/Users/ermay/Library/Fonts/Oswald-Bold.ttf')
oprop_r=fm.FontProperties(fname=opath_r)
oprop_sb=fm.FontProperties(fname=opath_sb)
oprop_b=fm.FontProperties(fname=opath_b)



def contour_presplot(x,y,data,center,label,savename):
    fig,axes=plt.subplots(figsize=(8.5,8))
    plt.gcf().subplots_adjust(bottom=0.12,top=0.95,left=0.12,right=0.76)

    plt.style.use('classic')
    fig.patch.set_facecolor('none')

    rcParams['axes.linewidth'] = 4.0

    rcParams['xtick.major.size'] = 15
    rcParams['xtick.major.width'] = 4
    rcParams['xtick.minor.size'] = 7
    rcParams['xtick.minor.width'] = 2

    rcParams['ytick.major.size'] = 15
    rcParams['ytick.major.width'] = 4
    rcParams['ytick.minor.size'] = 7
    rcParams['ytick.minor.width'] = 2
    

    minV=np.int((np.nanmin(data)))*1.0
    maxV=np.ceil(np.nanmax(data))
    
    if center==True:
        cbar=plt.cm.BrBG
        lim=np.nanmax([np.abs(minV),np.abs(maxV)])
        minV=-1.0*lim
        maxV=lim
    else:
        cbar=plt.cm.magma

    sig,lat=np.meshgrid(x,y)
    p=plt.contourf(sig,lat,data,vmin=minV,vmax=maxV,cmap=cbar)

    plt.xlabel('Latitude [degrees]',
               fontsize=25,fontproperties=oprop_b)
    plt.xticks(fontsize=20,fontproperties=oprop_r)
    plt.xlim(-90,90)

    plt.ylabel('Pressure [bar]',fontsize=25,fontproperties=oprop_b)
    plt.yticks(fontsize=20,fontproperties=oprop_r,
               rotation='vertical')
    plt.ylim(np.nanmax(y),np.nanmin(y))
    
    plt.figtext(0.78,0.12+(0.95-0.12)/2.0,label, 
            fontsize=25,fontproperties=oprop_b,rotation='vertical',
           ha='left',va='center')

    cbaxes = fig.add_axes([0.83, 0.12, 0.06, 0.95-0.12]) 
    cb = plt.colorbar(p,orientation='vertical',cax=cbaxes)
    cb.ax.tick_params(labelsize=20)
    cb.ax.set_xticklabels(cb.ax.get_xticks(),fontproperties=oprop_b)
    
    plt.show()
    plt.savefig(savename)
        