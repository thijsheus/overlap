#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:33:32 2019

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



###########################################################################

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
    
########################################################################
        
        
        
        
        
        
        
    line=np.arange(0,100,1)
    
    ###########################################
    ###### calculating overlap due to contrabution
    #fract= (5*np.pi) / (6*np.sqrt(3))  -1
    #fractal=  (5*np.pi) / (6*np.sqrt(3)) # -1
    #fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) -1
    
    shape=  (1/areaz) - 1
    shear=   (1/overlap_ratio) - (1/overlap)  
    #shear = abs(shear)
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
    
    

    FF = f_turb(ht,200);FF[ht<=100]=0;
    shear_sum = abs(shear_sum)
    shift_d = np.zeros(ht.size)
    shift_d[ht>50] = shift_distance
    
    shift_d = shift_d + shape*cp**0.5/2 ### draw the picture
    updraft = (shear_sum*ht)/shift_d
    updraft[updraft<0.4] = 0.4 ### based on neggers 2015
    FS = (shear_sum*ht)/(updraft*cp**0.5)
    #FS[updraft<0.1] = 0
    #FS[FS>max(1/overlap_ratio)] = np.nan
    #FS[np.argwhere(np.isnan(FS))] = max(FS)
    
    
    
    #shear = FS ; fractalf = FF
    
    r_total =  ( (shear) + (shape)  + (fractalf) ) +1
    
    
    #################################################################################################################
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  r_total -1
    vertical= 1/overlap_ratio -1
    
    
    bins=dz;
    #cmin=3
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
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        plt.xlim([0,13])
        plt.ylim([0,13])
        
    else:
        plt.clim(1,10**3)
        plt.xlim([0,10])
        plt.ylim([0,10])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_total_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_total_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_total_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_total_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_total_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    
    
    #plt.savefig('arm41400_total_emp1_r_061719.eps', dpi=300,bbox_inches='tight')
    """
    ###################################################################################
    """
    binwidth=0.1
    ###prob plot, histogram ,line plot
    plt.figure()
    #ybin, binEdges, patches =plt.hist(r_total*overlap_ratio,bins=10,density=True);#,log=True);  # may want to use one w/ log & one standard
    ybin, binEdges, patches =plt.hist(r_total*overlap_ratio,bins=np.arange(min(r_total*overlap_ratio), max(r_total*overlap_ratio) + binwidth, binwidth),density=True);
    plt.yticks(np.arange(0,10,2))
    
    plt.title('Accuracy of Model')
    plt.xlabel('Relative Effect')
    plt.ylabel('Probability Density')
    
    
    if m0==0:
        plt.savefig('rico82800_rel_effect_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_rel_effect_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_rel_effect_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_rel_effect_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_rel_effect_061919.eps', dpi=300,bbox_inches='tight')
    
    
    #plt.savefig('arm41400_rel_effect_061119.eps', dpi=300,bbox_inches='tight')
    """
    #################################################################################
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  shear
    vertical= 1/overlap_ratio -1
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=3)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=3)
    
    #plt.title('Effect of All Contributors');
    plt.title('Effect of Shear');
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
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        plt.xlim([0,10])
        plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        plt.xlim([0,8])
        plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_shear_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_shear_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_shear_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_shear_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_shear_r_061919.eps', dpi=300,bbox_inches='tight')
    
        
    
    
    #plt.savefig('arm41400_total_emp1_r_061719.eps', dpi=300,bbox_inches='tight')
    """
    #############################################################################################

    #################################################################################
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  shape
    vertical= 1/overlap_ratio -1
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=3)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=3)
    
    #plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    plt.title('Effect of Area Variability');
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
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        plt.xlim([0,10])
        plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        plt.xlim([0,8])
        plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_shape_r_061919.eps', dpi=300,bbox_inches='tight')
    
    
    
    
    #plt.savefig('arm41400_total_emp1_r_061719.eps', dpi=300,bbox_inches='tight')
    """
    #######################################################################################################
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  shear+shape
    vertical= 1/overlap_ratio -1
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=3)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=3)
    
    #plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    plt.title('Effect of Shear, Area Var.');
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
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        plt.xlim([0,10])
        plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        plt.xlim([0,8])
        plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_shear_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_shear_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_shear_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_shear_shape_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_shear_shape_r_061919.eps', dpi=300,bbox_inches='tight')
    """
    
    
    #######################################################################################
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  fractalf
    vertical= 1/overlap_ratio -1
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=3)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=3)
    
    #plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    #plt.title('Effect of Shear, Area Var.');
    plt.title('Effect of Turbulence');
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
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        plt.xlim([0,10])
        plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        plt.xlim([0,8])
        plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_turb_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_turb_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_turb_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_turb_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_turb_r_061919.eps', dpi=300,bbox_inches='tight')
    """
    ###################################################################################
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  1/overlap_ratio 
    vertical=  ht
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    #plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=3)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=3)
    
    #plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    #plt.title('Effect of Shear, Area Var.');
    #plt.title('Effect of Turbulence');
    plt.title('Cloud Height vs. Overlap');
    #plt.title('Overlap from Shear')
    #plt.title('Overlap from Turbulence')
    
    
    plt.xlabel('Actual Overlap')
    plt.ylabel('Cloud Height')
    #plt.xlabel('Numerical Overlap')
    #plt.ylabel('Overlap from Eqn')
    #plt.xlabel('Cloud Height')
    #plt.ylabel('Calculated Overlap')
    
    #colorbar=plt.colorbar(extend='both');
    colorbar=plt.colorbar()
    colorbar.set_label('count in bin')
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        #plt.xlim([0,10])
        #plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        #plt.xlim([0,8])
        #plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_ht_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_ht_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_ht_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_ht_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_ht_r_061919.eps', dpi=300,bbox_inches='tight')
    """

    ###################################################################################
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  1/overlap_ratio 
    vertical=  cp**0.5
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    #plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=3)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=3)
    
    #plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    #plt.title('Effect of Shear, Area Var.');
    #plt.title('Effect of Turbulence');
    plt.title('Cloud Width vs. Overlap');
    #plt.title('Overlap from Shear')
    #plt.title('Overlap from Turbulence')
    
    
    plt.xlabel('Actual Overlap')
    plt.ylabel('Cloud Width')
    #plt.xlabel('Numerical Overlap')
    #plt.ylabel('Overlap from Eqn')
    #plt.xlabel('Cloud Height')
    #plt.ylabel('Calculated Overlap')
    
    #colorbar=plt.colorbar(extend='both');
    colorbar=plt.colorbar()
    colorbar.set_label('count in bin')
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        #plt.xlim([0,10])
        #plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        #plt.xlim([0,8])
        #plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_width_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_width_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_width_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_width_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_width_r_061919.eps', dpi=300,bbox_inches='tight')
    """
    ##########################################################################################
    """
    #fract= (5*np.pi) / (6*np.sqrt(3))  -1
    fractal=  (5*np.pi) / (6*np.sqrt(3)) # -1
    fractalf_comb = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) -1
    
    
    
    #fractalf = np.zeros(ht.size)
    overlap_convex=overlap_convex[overlap_convex>0]
    turb= (1/overlap_ratio[ht>100]) - (1/overlap_convex)
    turb[turb<0]=0
    fractalf_comb[ht>100]=turb
    
    
    ## upper bound to turb
    lam =   (3*np.pi)/4  #np.pi * np.sqrt(3) / 2     #(15/8)*np.sqrt(3/8)*np.pi
    #C1=np.argwhere((1/overlap_ratio)>1.5*lam)
    keep=np.argwhere(fractalf_comb > lam)
    fractalf_comb[np.array(keep)] = lam
    
    
    #shear = FS ; fractalf = FF
    
    horizontal=  ht
    vertical=  fractalf_comb
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    #plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=5)
    plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=5)
    
    #plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    #plt.title('Effect of Shear, Area Var.');
    #plt.title('Effect of Turbulence');
    #plt.title('Overlap vs. Cloud Height');
    #plt.title('Overlap from Shear')
    plt.title('Overlap from Turbulence')
    
    
    #plt.xlabel('Calculated Overlap')
    #plt.ylabel('Actual Overlap')
    #plt.xlabel('Numerical Overlap')
    #plt.ylabel('Overlap from Equation')
    plt.xlabel('Cloud Height')
    plt.ylabel('Numerical Overlap')
    
    #colorbar=plt.colorbar(extend='both');
    colorbar=plt.colorbar()
    colorbar.set_label('count in bin')
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        #plt.xlim([0,10])
        #plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        #plt.xlim([0,8])
        #plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_turb_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_turb_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_turb_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_turb_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_turb_fit_r_061919.eps', dpi=300,bbox_inches='tight')
    """
    
    #########################################################################################################################
    
    """
    #shear = FS ; fractalf = FF
    
    horizontal=  shear
    vertical=  FS
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=5)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=5)
    
    #plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    #plt.title('Effect of Shear, Area Var.');
    #plt.title('Effect of Turbulence');
    #plt.title('Overlap vs. Cloud Height');
    plt.title('Overlap from Shear')
    #plt.title('Overlap from Turbulence')
    
    
    #plt.xlabel('Calculated Overlap')
    #plt.ylabel('Actual Overlap')
    plt.xlabel('Numerical Overlap')
    plt.ylabel('Overlap from Equation')
    #plt.xlabel('Cloud Height')
    #plt.ylabel('Numerical Overlap')
    
    #colorbar=plt.colorbar(extend='both');
    colorbar=plt.colorbar()
    colorbar.set_label('count in bin')
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        plt.xlim([0,6])
        plt.ylim([0,6])
        
    else:
        plt.clim(1,10**3)
        plt.xlim([0,4])
        plt.ylim([0,4])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    
    if m0==0:
        plt.savefig('rico82800_shear_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_shear_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_shear_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_shear_fit_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_shear_fit_r_061919.eps', dpi=300,bbox_inches='tight')
    
    
    """
    ##############################################################################################
    
    shear = FS ; fractalf = FF 
    r_total =  ( (shear) + (shape)  + (fractalf) ) +1
    
    horizontal=  r_total -1
    vertical=  1/overlap_ratio - 1
    
    
    bins=dz;
    #cmin=3
    ### cmin=3 for rico , concat bomex cmin=10 , concat bomex has many more cld's than 1 timestep of rico
    plt.figure();
    HIST2d=plt.hist2d(horizontal,vertical,bins=[1*dz,1*dz],cmin=cmin, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
    #plt.plot(horizontal,horizontal)
    #plt.plot(vertical,vertical)
    plt.plot(line,'k',linewidth=3)
    #plt.plot(ht,FF,'ko')
    
    #plt.plot(np.linspace(0,100,ht.size),np.zeros(ht.size),'k',linewidth=5)
    #plt.plot(np.linspace(100,max(ht)+25,ht.size),f_turb(np.linspace(100,max(ht)+25,ht.size),200),'k',linewidth=5)
    
    plt.title('Effect of All Contributors');
    #plt.title('Effect of Shear');
    #plt.title('Effect of Area Variability');
    #plt.title('Effect of Shear, Area Var.');
    #plt.title('Effect of Turbulence');
    #plt.title('Overlap vs. Cloud Height');
    #plt.title('Overlap from Shear')
    #plt.title('Overlap from Turbulence')
    
    
    plt.xlabel('Calculated Overlap')
    plt.ylabel('Actual Overlap')
    #plt.xlabel('Numerical Overlap')
    #plt.ylabel('Overlap from Equation')
    #plt.xlabel('Cloud Height')
    #plt.ylabel('Numerical Overlap')
    
    #colorbar=plt.colorbar(extend='both');
    colorbar=plt.colorbar()
    colorbar.set_label('count in bin')
    
    if m0==4: ### bomex
        plt.clim(1,10**4)
        plt.xlim([0,10])
        plt.ylim([0,10])
        
    else:
        plt.clim(1,10**3)
        plt.xlim([0,8])
        plt.ylim([0,8])
    
    
    
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
      
    #plt.loglog(horizontal,vertical,'o')
    
    """
    if m0==0:
        plt.savefig('rico82800_total_emp1_r_061919.eps', dpi=300,bbox_inches='tight')
        
    elif m0==1:
        plt.savefig('rico201600_total_emp1_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==2:
        plt.savefig('arm28800_total_emp1_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0==3:
        plt.savefig('arm41400_total_emp1_r_061919.eps', dpi=300,bbox_inches='tight')
        
        
    elif m0 ==4:
        plt.savefig('bomexAll_total_emp1_r_061919.eps', dpi=300,bbox_inches='tight')
    """


###############################################
end= time.time()
print('Run Time in Seconds:', end-start)




