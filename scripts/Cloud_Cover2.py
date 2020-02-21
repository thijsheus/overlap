#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:16:04 2019

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
#import os
#import sys
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
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

#######################################################

#### definitions
def f_turb(h,parm):
    return (3*np.pi /4)* (h/(parm+h))

def f_shear(aspect,parm):
    return parm*aspect**1

def f_area(aspect,parm):
    return parm*aspect**1


################################################
start=time.time()




conditional_height=0

file1_numb=-1

### accessing multiple datasets easily
Afilenames2 = sorted(glob.glob('/data/bomex/*.track.nc'))
#Bfilenames = Afilenames[5:]
Afilenames0 = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
Bfilenames = Afilenames0[22:23]
#Bfilenames = Afilenames[55:56]
Afilenames1 = sorted(glob.glob('/data/arm/*.track.nc'))
#Bfilenames = Afilenames[1:]
#Bfilenames = Afilenames[10:11]
#Bfilenames = Afilenames[17:18]
Cfilenames0 = [Afilenames0[22:23] , Afilenames0[55:56], Afilenames1[10:11], Afilenames1[17:18], Afilenames2[5:6]]

#Zfilenames = [Afilenames0,Afilenames1,Afilenames2]

m0=-1

for file1 in Cfilenames0:
    m0=m0+1
    file1=file1[0]
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
    
    if m0 == 0:
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
        cmin=3
        
    elif m0==1:
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
        cmin=3
        
    elif m0 ==2:
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
        cmin=3
        
    elif m0==3:
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
        cmin=3
        
    elif m0==4:
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
        cmin=5
        print('Bomex')
    
    
    """
    index1=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
    index1_size=index1[0].size
    ht=ht[index1[0]]   # taking the ht values according to indices above
    overlap_ratio=overlap_ratio[index1[0]];
    area_proj=area_proj[index1[0]]
    print('Clouds that satisfy the given condition: ',index1_size)
    print('conditional height is',conditional_height,'meters')
    """

    #fract= (5*np.pi) / (6*np.sqrt(3))  -1
    #fractal=  (5*np.pi) / (6*np.sqrt(3)) # -1
    #fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) -1
    
    shape=  (1/areaz) - 1
    shear=   (1/overlap_ratio) - (1/overlap)  
    shear[shear<0]=0
    
    fractalf = np.zeros(ht.size)
    overlap_convex=overlap_convex[overlap_convex>0]
    turb= (1/overlap_ratio[ht>100]) - (1/overlap_convex)
    turb[turb<0]=0
    fractalf[ht>100]=turb
    
    """
    ## upper bound to turb
    lam =   (3*np.pi)/4  #np.pi * np.sqrt(3) / 2     #(15/8)*np.sqrt(3/8)*np.pi
    #C1=np.argwhere((1/overlap_ratio)>1.5*lam)
    keep=np.argwhere(fractalf > lam)
    fractalf[np.array(keep)] = lam
    """
    
    r_total =  ( (shear) + (shape)  + (fractalf) ) +1
    
    
    #f_t=f_turb(ht,175); f_t[ht<=100]=0
    #f_s=f_shear(ht/l,0.43);f_s[ht/l<=1]=0
    #f_a=f_area(ht/l,0.2);f_a[ht/l<=1]=0
    
    #r_total = f_t + f_s + f_a +1
    
    cv_sorted=np.sort(cv)
    cv_arg=np.argsort(cv)
    r_sorted=np.zeros(cv.size)
    act_sorted=np.zeros(cv.size)
    count=-1
    for i in cv_arg:
        count= count+1
        r_sorted[count]=r_total[i]
        act_sorted[count]=1/overlap_ratio[i]
    
    
    ########
    #p_tot_trap=np.trapz(r_sorted,cv_sorted)
    #p_act_trap=np.trapz(act_sorted,cv_sorted)
    #print('Calculated:', p_tot_trap, 'Actual:',p_act_trap)
    #### total cld cover
    p_tot_num=sum(r_total*cv)
    p_act_num=sum(1/overlap_ratio * cv)
    print('Calculated:', p_tot_num, 'Actual:',p_act_num)

    

    #per_err= (abs(p_tot_trap - p_act_trap) / p_act_trap)*100
    #print(per_err)
    
    #perr=(p_tot_trap / p_act_trap) *100
    #print(perr)

    #per_err= (abs(p_tot_num - p_act_num) / p_act_num)*100
    #print(per_err)

    perr=(p_tot_num / p_act_num) *100
    print('percent of cld cover explained: ',perr)
    


    plt.figure()
    horizontal= (r_total)*cv
    vertical=  cp
    
    ## normalizing 
    max_cp = max(cp)
    horizontal = horizontal / max_cp
    vertical = vertical / max_cp
    
    plt.figure();
    HIST2dC=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=1,cmap='viridis',norm=colors.LogNorm()); #norm=colors.LogNorm()
    plt.plot(horizontal,horizontal,'k', linewidth=3)
    plt.title('Normalized Cloud Cover');

    
    plt.xlabel('Normalized Calculated Cloud Cover')
    plt.ylabel('Normalized Actual Cloud Cover')
    #plt.xlabel('Cloud Overlap')
    #plt.ylabel('Cloud Height')
    
    #colorbar=plt.colorbar(extend='both');
    colorbar=plt.colorbar()
    colorbar.set_label('count in bin')
    #plt.clim(0,200)
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        #plt.xlim([0,10])
        #plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        #plt.xlim([0,8])
        #plt.ylim([0,8])
        
    """
    if m0==0:
        plt.savefig('rico82800_cld_cover_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_cld_cover_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_cld_cover_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_cld_cover_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_cld_cover_061919.eps', dpi=300,bbox_inches='tight')
    """

    #plt.savefig('arm41400_cloud_cover_061119.eps', dpi=300, bbox_inches='tight')  
    
    

###############################################
end= time.time()
print('Run Time in Seconds:', end-start)


