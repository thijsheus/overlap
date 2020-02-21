#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 17:03:40 2019

@author: anthonys
"""



#from netCDF4 import Dataset
import numpy as np 
#import struct 
#import netCDF4
from netCDF4 import Dataset
#import netCDF4 as nc
#import collections
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
import matplotlib.colors as colors
#from scipy.io import netcdf
#import scipy as sp
import glob
import os
import sys
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import seaborn as sns
from scipy import stats
import time
#import pkgutil
#import collections
#from collections import Counter
from scipy.spatial import ConvexHull
import cProfile
from scipy.stats import spearmanr
from scipy.stats import kendalltau
#import overlap_calculation

plt.rcParams.update({'font.size': 16})

start=time.time()


################################################################################
# importing files

bomexd = Dataset("/data/bomex/bomex.default.0000000.nc","r")
bomexql = Dataset("/data/bomex/bomex.ql.0000000.nc","r")
bomexqlcore = Dataset("/data/bomex/bomex.qlcore.0000000.nc","r")
bomextrack18 = Dataset('/data/bomex/l.0001800.track.nc','r')
bomextrack36 = Dataset('/data/bomex/l.0003600.track.nc','r')
bomextrack54 = Dataset('/data/bomex/l.0005400.track.nc','r')
bomextrack72 = Dataset('/data/bomex/l.0007200.track.nc','r')
bomextrack90 = Dataset('/data/bomex/l.0009000.track.nc','r')
bomextrack108 = Dataset('/data/bomex/l.0010800.track.nc','r')
bomextrack126 = Dataset('/data/bomex/l.0012600.track.nc','r')
bomextrack144 = Dataset('/data/bomex/l.0014400.track.nc','r')
bomextrack162 = Dataset('/data/bomex/l.0016200.track.nc','r')
bomextrack180 = Dataset('/data/bomex/l.0018000.track.nc','r')
bomextrack198 = Dataset('/data/bomex/l.0019800.track.nc','r')
bomextrack216 = Dataset('/data/bomex/l.0021600.track.nc','r')

bomextrack342 = Dataset('/data/bomex/l.0034200.track.nc','r')
bomextrack360 = Dataset('/data/bomex/l.0036000.track.nc','r')


ricod = Dataset("/data/rico/rico.default.0000000.nc","r")
ricoql = Dataset("/data/rico/rico.ql.0000000.nc","r")
ricoqlcore = Dataset("/data/rico/rico.qlcore.0000000.nc","r")
ricotrack36 = Dataset('/data/rico/l.0003600.track.nc','r')
ricotrack72 = Dataset('/data/rico/l.0007200.track.nc','r')
ricotrack108 = Dataset('/data/rico/l.0010800.track.nc','r')
#ricotrack144 = Dataset('/data/rico/l.0014400.track.nc','r')
ricotrack180 = Dataset('/data/rico/l.0018000.track.nc','r')
ricotrack216 = Dataset('/data/rico/l.0021600.track.nc','r')
ricotrack252 = Dataset('/data/rico/l.0025200.track.nc','r')
ricotrack288 = Dataset('/data/rico/l.0028800.track.nc','r')
ricotrack324 = Dataset('/data/rico/l.0032400.track.nc','r')
ricotrack360 = Dataset('/data/rico/l.0036000.track.nc','r')
ricotrack396 = Dataset('/data/rico/l.0039600.track.nc','r')

ricotrack612 = Dataset('/data/rico/l.0061200.track.nc','r')
ricotrack828 = Dataset('/data/rico/l.0082800.track.nc','r')
ricotrack900 = Dataset('/data/rico/l.0090000.track.nc','r')
ricotrack1008 = Dataset('/data/rico/l.0100800.track.nc','r')
ricotrack1116 = Dataset('/data/rico/l.0111600.track.nc','r')
ricotrack1224 = Dataset('/data/rico/l.0122400.track.nc','r')
ricotrack1332 = Dataset('/data/rico/l.0133200.track.nc','r')
ricotrack1440 = Dataset('/data/rico/l.0144000.track.nc','r')
ricotrack1548 = Dataset('/data/rico/l.0154800.track.nc','r')
ricotrack1656 = Dataset('/data/rico/l.0165600.track.nc','r')
ricotrack1764 = Dataset('/data/rico/l.0176400.track.nc','r')
ricotrack1872 = Dataset('/data/rico/l.0187200.track.nc','r')
ricotrack1980 = Dataset('/data/rico/l.0198000.track.nc','r')


ricotrack2016 = Dataset('/data/rico/l.0201600.track.nc','r')
ricotrack2052 = Dataset('/data/rico/l.0205200.track.nc','r')
ricotrack2088 = Dataset('/data/rico/l.0208800.track.nc','r')
ricotrack2124 = Dataset('/data/rico/l.0212400.track.nc','r')
ricotrack2160 = Dataset('/data/rico/l.0216000.track.nc','r')



armd = Dataset("/data/arm/arm.default.0000000.nc","r")
armql = Dataset("/data/arm/arm.ql.0000000.nc","r")
armqlcore = Dataset("/data/arm/arm.qlcore.0000000.nc","r")
armtrack108 = Dataset('/data/arm/l.0010800.track.nc','r')
armtrack126 = Dataset('/data/arm/l.0012600.track.nc','r')
armtrack144 = Dataset('/data/arm/l.0014400.track.nc','r')
armtrack162 = Dataset('/data/arm/l.0016200.track.nc','r')
armtrack180 = Dataset('/data/arm/l.0018000.track.nc','r')
armtrack198 = Dataset('/data/arm/l.0019800.track.nc','r')
armtrack216 = Dataset('/data/arm/l.0021600.track.nc','r')
armtrack234 = Dataset('/data/arm/l.0023400.track.nc','r')
armtrack252 = Dataset('/data/arm/l.0025200.track.nc','r')
armtrack270 = Dataset('/data/arm/l.0027000.track.nc','r')
armtrack288 = Dataset('/data/arm/l.0028800.track.nc','r')

armtrack504 = Dataset('/data/arm/l.0050400.track.nc','r')
armtrack522 = Dataset('/data/arm/l.0052200.track.nc','r')



lassotrack306 = Dataset('/data/lasso/sims/20160611_micro/l.0030600.track.nc','r')

"""
filenames=[bomexd, ricod, armd]

bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72, bomextrack90, bomextrack108, bomextrack126, bomextrack144, bomextrack162, bomextrack180, bomextrack198, bomextrack216]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144, ricotrack180, ricotrack216, ricotrack252, ricotrack288, ricotrack324, ricotrack360, ricotrack396]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162, armtrack180, armtrack198, armtrack216, armtrack234, armtrack252, armtrack270, armtrack288]
"""
###########################################################################
# def 
"""
def overlap(s,h,l):
    #return l / (l +s*h)
    return 1 - (s*h / (l + s*h))

def overlap1(s,h,l,w):
    alpha= 1 / (1+w) + 1
    return 1 - (alpha*0.5*s*h / (l + s*h))
"""
#def overlap2(s,h,l):
 #   return 1 - (2*s*h / (np.pi*l + 2*s*h))

#### definitions
def f_turb(h,parm):
    return (3*np.pi /4)* (h/(parm+h))

def f_shear(aspect,parm):
    return parm*aspect**1

def f_area(aspect,parm):
    return parm*aspect**1

def f_dsdz(aspect,shear):
    return shear*aspect

####################################
    
# script
filenames=[ricod]
bomexfilenames=[bomextrack360]
ricofilenames=[ricotrack828]
#armfilenames=[armtrack522]
lassofilenames=[lassotrack306]

### rico828
"""
height=np.load('height_no_shear_rico828.npy')
volume=np.load('volume_no_shear_rico828.npy')
projar=np.load('projarea_no_shear_rico828.npy')
overlap=np.load('overlap_no_shear_rico828.npy')
areaz=np.load('area_z_ratio_rico828.npy')
areazg=np.load('rico_area_gmean82800.npy')
overlap_convex=np.load('overlap_convex_rico82800.npy')
#shearh_vel = np.load('rico82800_shearh_vel.npy')
#sheara_vel = np.load('rico82800_sheara_vel.npy')
#wavg = np.load('rico82800_wavg.npy')
#shearTB_vel = np.load('rico82800_shearTB_vel.npy')
wavg1 = np.load('wavg1_rico82800.npy')
shearTB = np.load('shearTB_rico82800.npy')
shear_sum = np.load('shear_sum_rico82800.npy')
wz_max = np.load('wz_max_rico82800.npy')
wz_cb = np.load('wz_cb_rico82800.npy')
shift_distance = np.load('shift_distancerico82800.npy')
"""
"""
### rico2016
height=np.load('height_no_shear_rico2016.npy')
volume=np.load('volume_no_shear_rico2016.npy')
projar=np.load('projarea_no_shear_rico2016.npy')
overlap=np.load('overlap_no_shear_rico2016.npy')
areaz=np.load('area_z_ratio_rico2016.npy')
model=np.load('model_rico2016.npy')
model_inv=np.load('model_inv_rico2016.npy') 
model_opt=np.load('model_opt_rico2016.npy')
model_opt2=np.load('model_opt2_rico2016.npy')
overlap_convex=np.load('overlap_convex_rico201600.npy')
wavg1 = np.load('wavg1_rico201600.npy')
shearTB = np.load('shearTB_rico201600.npy')
shear_sum = np.load('shear_sum_rico201600.npy')
shift_distance = np.load('shift_distancerico201600.npy')
"""
"""
### bomex 342
height342=np.load('height_no_shear_bomex34200.npy')
volume342=np.load('volume_no_shear_bomex34200.npy')
projar342=np.load('projarea_no_shear_bomex34200.npy')
overlap342=np.load('overlap_no_shear_bomex34200.npy')
areaz342=np.load('area_z_ratio_bomex34200.npy')

### bomex 360
height=np.load('height_no_shear_bomex36000.npy')
volume=np.load('volume_no_shear_bomex36000.npy')
projar=np.load('projarea_no_shear_bomex36000.npy')
overlap=np.load('overlap_no_shear_bomex36000.npy')
areaz=np.load('area_z_ratio_bomex36000.npy')

### lasso 306
height=np.load('height_no_shear_lasso306.npy')
volume=np.load('volume_no_shear_lasso306.npy')
projar=np.load('projarea_no_shear_lasso306.npy')
overlap=np.load('overlap_no_shear_lasso306.npy')
areaz=np.load('area_z_ratio_lasso306.npy')
"""
"""
### arm28800

height=np.load('height_no_shear_arm28800.npy')
volume=np.load('volume_no_shear_arm28800.npy')
projar=np.load('projarea_no_shear_arm28800.npy')
overlap=np.load('overlap_no_shear_arm28800.npy')
areaz=np.load('area_z_ratio_arm28800.npy')
overlap_convex=np.load('overlap_convex_arm28800.npy')
#shearh_vel = np.load('rico82800_shearh_vel.npy')
#sheara_vel = np.load('rico82800_sheara_vel.npy')
#wavg = np.load('rico82800_wavg.npy')
#shearTB_vel = np.load('rico82800_shearTB_vel.npy')
wavg1 = np.load('wavg1_arm28800.npy')
shearTB = np.load('shearTB_vel_arm28800.npy')
shear_sum = np.load('shear_sum_arm28800.npy')
shift_distance = np.load('shift_distancearm28800.npy')
"""

### arm41400

height=np.load('height_no_shear_arm41400.npy')
volume=np.load('volume_no_shear_arm41400.npy')
projar=np.load('projarea_no_shear_arm41400.npy')
overlap=np.load('overlap_no_shear_arm41400.npy')
areaz=np.load('area_z_ratio_arm41400.npy')
overlap_convex=np.load('overlap_convex_arm41400.npy')
#shearh_vel = np.load('rico82800_shearh_vel.npy')
#sheara_vel = np.load('rico82800_sheara_vel.npy')
#wavg = np.load('rico82800_wavg.npy')
#shearTB_vel = np.load('rico82800_shearTB_vel.npy')
wavg1 = np.load('wavg1_arm41400.npy')
shearTB = np.load('shearTB_vel_arm41400.npy')
shear_sum = np.load('shear_sum_arm41400.npy')
shift_distance = np.load('shift_distancearm41400.npy')

### Bomex_All is Below



conditional_height=0

file1_numb=-1

### accessing multiple datasets easily
#Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
#Bfilenames = Afilenames[5:]
#Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
#Bfilenames = Afilenames[22:23]
#Bfilenames = Afilenames[55:56]
Afilenames = sorted(glob.glob('/data/arm/*.track.nc'))
#Bfilenames = Afilenames[1:]
#Bfilenames = Afilenames[10:11]
Bfilenames = Afilenames[17:18]

for file1 in Bfilenames:
    print(file1)

    data = Dataset(file1,'r')

    file1_numb = file1_numb+1
    
    ht=data.variables['ht'][:]
    cb=data.variables['cb'][:]
    ct=data.variables['ct'][:]
    cv=data.variables['cv'][:]
    cp=data.variables['cp'][:]
    overlap_ratio=data.variables['chr'][:]
    area_proj=data.variables['area_proj'][:]
    nrcloud=data.variables['nrcloud'][:,:,:]
    cfrac=data.variables['cfrac'][:]
    zt=data.variables['z'][:]
    xt=data.variables['x'][:]
    yt=data.variables['y'][:]
    nr=data.variables['nr'][:]
    cld_mask=data.variables['cld_mask'][:,:,:]
    
    
    nrcloudarray = np.ma.getdata(nrcloud) # unmask array

    dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
    gridarea=dx*dy
    gridvol=dx*dy*dz
    nx=xt.size;ny=yt.size;nz=zt.size;

    l=cv**0.5
########################################################################
"""
### Bomex_All
overlap_ratio=np.load('Bomex_overlap_ratio.npy')            
height=np.load('Bomex_ht.npy')
overlap=np.load('Bomex_shear.npy')
areaz=np.load('Bomex_areaz.npy')
ht=height
cv=np.load('Bomex_cv.npy')
cp=np.load('Bomex_cp.npy')
overlap_convex=np.load('Bomex_hull.npy')
#shearh_vel = np.load('Bomex_shearh_vel.npy')
#sheara_vel = np.load('Bomex_sheara_vel.npy')
#wavg = np.load('Bomex_wavg.npy')
shearTB = np.load('Bomex_shearTB.npy')
shear_sum = np.load('Bomex_shear_sum.npy')
wz_max = np.load('Bomex_wz_max.npy')
wz_cb = np.load('Bomex_wz_cb.npy')
shift_distance = np.load('Bomex_shift_distance.npy')
print('Bomex')
"""


"""

bins=dz;
plt.figure();
plt.hist2d(overlap_c,height_c,bins=bins,cmin=0.5);
plt.title('Actual Overlap vs. Calculated Overlap with No Shear');
colorbar=plt.colorbar();
colorbar.set_label('count in bin')
"""


## quick calculation
#cross_area= np.extract(area_z_ratio <1, area_z_ratio)




line=np.arange(0,100,1)
#areaz=areazg

###########################################
###### calculating overlap due to contrabution
fract= (5*np.pi) / (6*np.sqrt(3))  -1
fractal=  (5*np.pi) / (6*np.sqrt(3)) # -1
fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) -1

shape=  (1/areaz) - 1
shear=   (1/overlap_ratio) - (1/overlap)  
#shear = abs(shear)
shear[shear<0]=0

overlap_convex=overlap_convex[overlap_convex>0]
turb= (1/overlap_ratio[ht>100]) - (1/overlap_convex)
turb[turb<0]=0
fractalf[ht>100]=turb


## upper bound to turb
lam =   (3*np.pi)/4  #np.pi * np.sqrt(3) / 2     #(15/8)*np.sqrt(3/8)*np.pi
#C1=np.argwhere((1/overlap_ratio)>1.5*lam)
keep=np.argwhere(fractalf > lam)
fractalf[np.array(keep)] = lam

##fractalf = np.minimum(1.5*lam*np.ones(overlap_ratio.size),1/overlap_ratio) -1 ### trying only 3d Koch overpredicts

#fractalf = np.minimum(((3*np.pi) / 4)*np.ones(overlap_ratio.size),1/overlap_ratio) -1

#fractalf = np.minimum( (3*np.pi /4)*np.ones(overlap_ratio.size), 1/overlap_ratio ) 
#fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) 

"""
fractal_ones=fractal*np.ones(overlap_ratio.size)
r_fractal= np.minimum(1/fractal_ones,1/overlap_ratio)
r_total =  ( (1/(overlap)) + (1/areaz) + (r_fractal) ) -3  # correction: need to subtract 1 from each 

plt.figure()
plt.plot(r_fractal,1/overlap_ratio,'o')
plt.plot(np.arange(0,3,1),'k',linewidth=3)
plt.xlim(left=0)
plt.ylim(bottom=0)
"""

#shape = 1 / np.exp( -ht/200 ) - 1

#r_total =  ( (1/(overlap)) + (1/areaz) + (fractal) )   
#r_total =  ( (shear) + (shape)  + (fractalf) ) +1  

#updraft_cb=np.mean(wz_cb)*np.zeros(ht.size)
#updraft_cb[ht>50] = wz_cb
#updraft_max=np.mean(wz_cb)*np.zeros(ht.size)
#updraft_max[ht>50] = wz_cb


FF = f_turb(ht,200);FF[ht<=100]=0;
shear_sum = abs(shear_sum)
shift_d = np.zeros(ht.size)
shift_d[ht>50] = shift_distance

shift_d = shift_d + shape*cp**0.5/2 ### draw the picture
updraft = (shear_sum*ht)/shift_d
FS = (shear_sum*ht)/(updraft*cp**0.5)
FS[updraft<0.1] = 0

#shear = FS ; fractalf = FF

r_total =  ( (shear) + (shape)  + (fractalf) ) +1  

horizontal=  r_total
vertical= 1/overlap_ratio


bins=dz;
cmin=3
### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
plt.figure();
HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
#plt.plot(horizontal,horizontal)
#plt.plot(vertical,vertical)
plt.plot(line,'k',linewidth=3)
#plt.plot(ht,FF,'ko')

#plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=3)
#plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=3)

plt.title('Effect of All Contributors');
#plt.title('Effect of Shear');
#plt.title('Effect of Area Variability');
#plt.title('Effect of Shear, Area Var.');
#plt.title('Effect of Turbulence');
#plt.title('Relating Overlap and Height');
#plt.title('Overlap from Shear')
#plt.title('Overlap from Turbulence')


plt.xlabel('Calculated Overlap')
plt.ylabel('Actual Overlap')
#plt.xlabel('Numerical Overlap')
#plt.ylabel('Overlap from Eqn')
#plt.xlabel('Cloud Height')
#plt.ylabel('Calculated Overlap')

#colorbar=plt.colorbar(extend='both');
colorbar=plt.colorbar()
colorbar.set_label('count in bin')
#plt.clim(0,200)

#plt.xlim([0,1])
#plt.ylim([0,1])
#plt.xlim(left=0)
#plt.ylim(bottom=0)
plt.xlim([0,8])
plt.ylim([0,8])  #use cmin=5
#plt.loglog(horizontal,vertical,'o')




#plt.savefig('arm41400_total_emp1_r_061719.eps', dpi=300,bbox_inches='tight')

###################################################################################
binwidth=0.1
###prob plot, histogram ,line plot
plt.figure()
#ybin, binEdges, patches =plt.hist(r_total*overlap_ratio,bins=10,density=True);#,log=True);  # may want to use one w/ log & one standard
ybin, binEdges, patches =plt.hist(r_total*overlap_ratio,bins=np.arange(min(r_total*overlap_ratio), max(r_total*overlap_ratio) + binwidth, binwidth),density=True);
plt.yticks(np.arange(0,10,2))

"""
#plt.figure()
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
plt.plot(bincenters,ybin,'k-',linewidth=2)
"""
plt.title('Accuracy of Model')
plt.xlabel('Relative Effect')
plt.ylabel('Probability Density')


#plt.savefig('arm41400_rel_effect_061119.eps', dpi=300,bbox_inches='tight')

#############################
"""
plt.figure()
plt.hist2d(ht,vertical-horizontal,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis', norm=colors.LogNorm());

plt.plot(ht,np.zeros(ht.size),'k',linewidth=3)


plt.xlabel('Cloud Height')
plt.ylabel('Residual')

colorbar=plt.colorbar()
colorbar.set_label('count in bin')
"""
################################################

"""
plt.figure()
# Density Plot and Histogram 
sns.distplot(r_total*overlap_ratio, hist=True, kde=False, 
             bins=int(10), color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 2})
"""

#plt.figure()
#sns.distplot(r_total*overlap_ratio)



##################################################################################
"""
binwidth=500000
###prob plot, histogram ,line plot
plt.figure()
horiz = (1/overlap) * cv
#ybin, binEdges, patches =plt.hist(r_total*overlap_ratio,bins=10,density=True);#,log=True);  # may want to use one w/ log & one standard
ybin, binEdges, patches =plt.hist(horiz,bins=np.arange(min(horiz), max(horiz) + binwidth, binwidth),density=False);
#plt.yticks(np.arange(0,10,2))


plt.title('Calculated Cloud Cover')
#plt.xlabel('Relative Effect')
#plt.ylabel('Probability Density')


cld_cover= sum(ybin*binwidth)
print(cld_cover)
"""
########################################################################################
"""
plt.figure()

plt.plot((1/overlap)-1,(1/overlap_ratio),'ob',label='Shear')
#plt.plot(line)
plt.plot((1/overlap)+(1/areaz)-2,(1/overlap_ratio),color='purple',marker='o',linestyle='none',label='Shear and Shape')
plt.plot((1/overlap)+(1/areaz)+(1/fractal)-3,(1/overlap_ratio),color='pink',marker='o',linestyle='None',label='Shear, Shape, and Fractal')
plt.plot(np.arange(0,10,1),'k',linewidth=3)
plt.xlabel('Inverse Calculated Overlap')
plt.ylabel('Inverse Actual Overlap')
plt.title('Inverse Overlap due to Contributions')
plt.legend(loc='lower right', bbox_to_anchor=(1.10, 0.0),fontsize='small',ncol=1)
#plt.savefig('rico2016_contributions_041619.eps', dpi=300, bbox_inches='tight')
"""
"""
plt.figure()
plt.subplot(2,2,1)
plt.hist2d(1/overlap -1,1/overlap_ratio,bins=bins,cmin=2,cmap='viridis');
plt.plot(line,'k',linewidth=3)
plt.subplot(2,2,2)
plt.hist2d( 1/areaz -1,1/overlap_ratio,bins=bins,cmin=2,cmap='viridis');
plt.plot(line,'k',linewidth=3)

plt.subplot(2,2,3)
plt.hist2d(1/np.ones(overlap_ratio.size)*fractal -1,1/overlap_ratio,bins=bins,cmin=2,cmap='viridis');
plt.plot(line,'k',linewidth=3)

plt.subplot(2,2,4)
plt.hist2d(1/overlap + 1/areaz + 1/fractal -3,1/overlap_ratio,bins=bins,cmin=2,cmap='viridis');
plt.plot(line,'k',linewidth=3)
"""
###################################################
##subplots 
"""
plt.figure()
fig, ax = plt.subplots(nrows=1, ncols=3, sharex=True, sharey=True, figsize=(9, 4));None

plt.subplot(1,3,1)
plt.hist2d(1/overlap -1,1/overlap_ratio,bins=bins,cmin=2,cmap='viridis');
plt.plot(line,'k',linewidth=3)
plt.clim(0,200)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.title('Shear')

plt.subplot(1,3,2)
plt.hist2d( 1/overlap + 1/areaz -2,1/overlap_ratio,bins=bins,cmin=2,cmap='viridis');
plt.plot(line,'k',linewidth=3)
plt.clim(0,200)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.title('Shear,Area')

plt.subplot(1,3,3)
plt.hist2d(1/overlap + 1/areaz + 1/fractal -3,1/overlap_ratio,bins=bins,cmin=2,cmap='viridis');
plt.plot(line,'k',linewidth=3)
plt.clim(0,200)
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.title('All')

#plt.suptitle('Inverse Overlap');
colorbar=plt.colorbar(extend='both');
colorbar.set_label('count in bin')
plt.clim(0,200)
fig.text(0.5, 0.004, 'Inverse Calculated Overlap', ha='center')
fig.text(0.06, 0.5, 'Inverse Actual Overlap', va='center', rotation='vertical')

#plt.xlabel('Inverse Calculated Overlap')
#plt.ylabel('Inverse Actual Overlap')
#plt.xlim([0,1])
#plt.ylim([0,1])


#plt.savefig('rico2016_subplot_invoverlap_042219.eps', dpi=300, bbox_inches='tight')
"""
############################
"""
plt.figure()

plt.figure();
plt.hist2d(1/overlap -1, 1/overlap_ratio,bins=bins,cmin=2,cmap='pink',alpha=0.3); #norm=colors.LogNorm(), normed: normalizes data
plt.hist2d(1/overlap + 1/areaz -2, 1/overlap_ratio,bins=bins,cmin=2,cmap='bone',alpha=0.3);
plt.hist2d(1/overlap -1 + 1/areaz + 1/fractal -3, 1/overlap_ratio,bins=bins,cmin=2,cmap='Blues',alpha=0.3);
#plt.plot(horizontal,horizontal)
#plt.plot(vertical,vertical)
plt.plot(line,'k',linewidth=3)
#plt.title('Actual Overlap vs. Area Ratio');
plt.title('Inverse Overlap');
colorbar=plt.colorbar(extend='both');
colorbar.set_label('count in bin')
plt.clim(0,200)
plt.xlabel('Inverse Total Calculated Overlap')
plt.ylabel('Inverse Actual Overlap')
"""

"""
plt.figure()
plt.scatter(1/overlap -1, 1/overlap_ratio,cmap='pink',alpha=0.3)
plt.scatter(1/overlap + 1/areaz -2, 1/overlap_ratio,cmap='bone',alpha=0.3);
plt.scatter(1/overlap -1 + 1/areaz + 1/fractal -3, 1/overlap_ratio,cmap='Blues',alpha=0.3);
plt.plot(np.arange(0,10,1),'k',linewidth=3)
"""











###############################################
end= time.time()
print('Run Time in Seconds:', end-start)











































