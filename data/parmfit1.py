#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 11:32:06 2019

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
#from scipy.io import netcdf
#import scipy as sp
import glob
#import os
#import sys
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.colors as colors
import time
#import pkgutil
#import collections
#from collections import Counter
from scipy.spatial import ConvexHull
#import cProfile
from numpy import percentile
#import overlap_calculation
from scipy.stats import spearmanr
#from scipy.stats import hmean
from scipy.stats import gmean



################################################
start=time.time()


def shearpar(X,gs):
    h,l = X
    return gs * (h/l)


def shearpar2(X,gs):
    h,l = X
    return gs * (h/l)**2


def turbpar(h,g):
    return g*h

def multi1(X,parm1,parm2,parm3):
    h,l,V = X
    return parm1*h+parm2*l+parm3*V


def f_turb(h,parm):
    return (3*np.pi /4)* (h/(parm+h))


####################################################

### rico828
height=np.load('height_no_shear_rico828.npy')
volume=np.load('volume_no_shear_rico828.npy')
projar=np.load('projarea_no_shear_rico828.npy')
overlap=np.load('overlap_no_shear_rico828.npy')
areaz=np.load('area_z_ratio_rico828.npy')
areazg=np.load('rico_area_gmean82800.npy')
overlap_convex=np.load('overlap_convex_rico82800.npy')

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
"""





conditional_height=0

file1_numb=-1

### accessing multiple datasets easily
#Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
#Bfilenames = Afilenames[5:]
#Bfilenames = Afilenames[5:6]
Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
Bfilenames = Afilenames[22:23]
#Bfilenames = Afilenames[55:56]

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
    
    """
    ### Bomex 
    overlap_ratio=np.load('Bomex_overlap_ratio.npy')            
    height=np.load('Bomex_ht.npy')
    overlap=np.load('Bomex_shear.npy')
    areaz=np.load('Bomex_areaz.npy')
    ht=height
    cv=np.load('Bomex_cv.npy')
    areazg=np.load('Bomex_area_gmean.npy')
    overlap_convex=np.load('Bomex_hull.npy')
    print('Bomex')
    """
    
    """
    index1=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
    index1_size=index1[0].size
    ht=ht[index1[0]]   # taking the ht values according to indices above
    overlap_ratio=overlap_ratio[index1[0]];
    area_proj=area_proj[index1[0]]
    print('Clouds that satisfy the given condition: ',index1_size)
    print('conditional height is',conditional_height,'meters')
    """
    
    
    
    fract= (5*np.pi) / (6*np.sqrt(3))  -1
    fractal=  (5*np.pi) / (6*np.sqrt(3)) # -1
    fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) -1
    
    shape=  (1/areaz) - 1
    shear=   (1/overlap_ratio) - (1/overlap)  
    shear[shear<0]=0
    
    turb= (1/overlap_ratio[ht>100]) - (1/overlap_convex)
    turb[turb<0]=0
    fractalf[ht>100]=turb
    """
    ### upper bound to turb
    lam =  np.pi * np.sqrt(3) / 2     #(15/8)*np.sqrt(3/8)*np.pi
    keep=np.argwhere(fractalf >1.5*lam-1)
    fractalf[np.array(keep)]=1.5*lam-1
    """
    
    ### upper bound to turb
    lam =   (3*np.pi)/4  #np.pi * np.sqrt(3) / 2     #(15/8)*np.sqrt(3/8)*np.pi
    #C1=np.argwhere((1/overlap_ratio)>1.5*lam)
    keep=np.argwhere(fractalf > lam)
    fractalf[np.array(keep)] = lam
    
    r_total =  ( (shear) + (shape)  + (fractalf) ) +1
    
    l=cv**0.5
    
    """
    cmin=1
    l=cv**0.5
    V=volume
    X1=(ht,l,V)
    opt1, cov1 = curve_fit(multi1, X1 , shear)
    model_opt1=multi1(X1,*opt1)
    
    plt.figure()
    plt.hist2d(shear, model_opt1, bins=dz , cmin=cmin)
    plt.plot(np.arange(0,100,1),'k',linewidth=3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    #plt.xlim([0,20])
    plt.xlabel('Calculated Overlap')
    plt.ylabel('Actual Overlap')
    colorbar=plt.colorbar(extend='both');
    colorbar.set_label('count in bin')
    plt.clim(0,200)
    
    
    cmin=1
    l=cv**0.5
    X1=(ht,l)
    opt1, cov1 = curve_fit(shearpar, X1 , shear,p0=2)
    model_opt1=shearpar(X1,*opt1)
    
    plt.figure()
    plt.hist2d(ht/l, model_opt1, bins=dz , cmin=cmin)
    plt.plot(np.arange(0,100,1),'k',linewidth=3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    #plt.xlim([0,20])
    plt.xlabel('Calculated Overlap')
    plt.ylabel('Actual Overlap')
    colorbar=plt.colorbar(extend='both');
    colorbar.set_label('count in bin')
    plt.clim(0,200)
    
    
    
    
    ht_sort=np.sort(ht)
    ht_arg=np.argsort(ht)
    fractalf_sort=np.zeros(fractalf.size)
    count=-1
    for i in ht_arg:
        count= count+1
        fractalf_sort[count]=fractalf[i]
    
    cmin=1
    l=cv**0.5
    X2=ht
    opt2, cov2 = curve_fit(turbpar, X2 , fractalf)
    model_opt2=turbpar(X2,*opt2)
    
    plt.figure()
    plt.hist2d(ht,model_opt2, bins=dz , cmin=cmin)
    plt.plot(np.arange(0,100,1),'k',linewidth=3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    #plt.xlim([0,20])
    plt.xlabel('Calculated Overlap')
    plt.ylabel('Actual Overlap')
    colorbar=plt.colorbar(extend='both');
    colorbar.set_label('count in bin')
    plt.clim(0,200)
    """
    
    """
    ### curve fit binned data
    opt_t, cov_t = curve_fit(f_turb , ht, fractalf)

    c=opt_t
    model_turb=f_turb(ht,c)
    
    
    bins=1*dz;
    
    plt.figure();
    plt.hist2d(ht,model_turb,bins=[1*dz,1*dz],cmin=1, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    
    
    plt.title('Effect of Turbulence');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    #plt.title('Effect of Shear, Area Var.');
    #plt.title('Effect of Turbulence');
    #plt.title('Relating Overlap and Height');
    
    plt.xlabel('ht')
    plt.ylabel('Cal Overlap')
    #plt.xlabel('Cloud Overlap')
    #plt.ylabel('Cloud Height')
    
    #colorbar=plt.colorbar(extend='both');
    colorbar=plt.colorbar()
    colorbar.set_label('count in bin')
    """
    
    plt.figure()
    plt.plot(ht,fractalf,'o')
    z= np.polyfit(ht, fractalf,3)
    p=np.poly1d(z)
    plt.plot(ht,p(ht),'o')
    
    
    plt.figure()
    plt.plot(fractalf,p(ht),'o')
    plt.plot(fractalf,fractalf)
    
    
    #plt.figure()
    ##plt.plot(fractalf-p(ht),'o') # residual plot
    
    #plt.hist(fractalf-p(ht))
    
    plt.figure()
    plt.hist2d(ht,fractalf - p(ht), bins=25,cmin=1) # res plt
    plt.clim(0,200)
    plt.plot(ht,np.zeros(ht.size),'k',linewidth=3)
    
    z_t = z;p_t=p
    
    plt.figure()
    plt.plot(ht/l,shear,'o')
    z= np.polyfit(ht/l, shear,3)
    p=np.poly1d(z)
    plt.plot(ht/l,p(ht/l),'o')
    
    
    plt.figure()
    plt.plot(shear,p(ht/l),'o')
    plt.plot(shear,shear) 
    
    plt.figure()
    plt.hist2d(ht/l,shear - p(ht/l), bins=25,cmin=1) # res plt
    plt.clim(0,200)
    plt.plot(ht/l,np.zeros(ht.size),'k',linewidth=3)
    z_s = z;p_s=p
    
    
    plt.figure()
    plt.plot(ht,shape,'o')
    z= np.polyfit(ht, shape,3)
    p=np.poly1d(z)
    plt.plot(ht,p(ht),'o')
    
    
    plt.figure()
    plt.plot(shape,p(ht),'o')
    plt.plot(shape,shape)
    
    plt.figure()
    plt.hist2d(ht,shape - p(ht), bins=25,cmin=1) # res plt
    plt.clim(0,200)
    plt.plot(ht,np.zeros(ht.size),'k',linewidth=3)
    z_a = z;p_a=p
    
    
    
    p_c = p_s(ht/l) + p_t(ht) + p_a(ht) +1
    
    
    
    plt.figure()
    plt.plot(p_c,1/overlap_ratio,'o')
    plt.plot(p_c,p_c) 
    
    plt.figure()
    plt.hist2d(p_c,1/overlap_ratio, bins=dz , cmin=3)
    plt.plot(np.arange(0,100,1),'k',linewidth=3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    #plt.xlim([0,20])
    plt.xlabel('Calculated Overlap')
    plt.ylabel('Actual Overlap')
    colorbar=plt.colorbar(extend='both');
    colorbar.set_label('count in bin')
    plt.clim(0,200)



###############################################
end= time.time()
print('Run Time in Seconds:', end-start)