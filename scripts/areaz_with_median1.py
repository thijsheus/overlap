#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:24:46 2019

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


conditional_height=0

file1_numb=-1

### accessing multiple datasets easily
begin=5
### accessing multiple datasets easily
Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
Bfilenames = Afilenames[begin:]
#Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
#Bfilenames = Afilenames[22:23]

file1numb=begin-1


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



    index1=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
    index1_size=index1[0].size
    ht=ht[index1[0]]   # taking the ht values according to indices above
    overlap_ratio=overlap_ratio[index1[0]];
    area_proj=area_proj[index1[0]]
    print('Clouds that satisfy the given condition: ',index1_size)
    print('conditional height is',conditional_height,'meters')


    cloud_numb=index1[0] +1 # index is off by 1 as array starts with zero
    
    alpha=np.arange(1,0,-1) # just 1
    overlap_changed_matrix=np.zeros((cloud_numb.size,alpha.size))
    area_z_ratio=np.zeros(cloud_numb.size)
    m=-1
    for e1 in cloud_numb:
        m=m+1
        location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1
        base=np.amin(location[0]);top=np.amax(location[0]);
        layers=(top - base + 1)   # layers of cloud along z 
        cross_area_z=np.zeros(int(layers))
        

        for z in np.arange(base,top+1,1):
            findz=np.where(location[0] == z)
            cross_area_z[z-base]=dx*dy*len(findz[0])  
            
            
        area_z_ratio[m]= gmean(cross_area_z) / np.amax(cross_area_z)
        
    step=1800 # 30*60
    time_array= np.arange(10800,36000+step,step)
    npyfilespath ="/home/anthonys/Documents/bomex_area_gmean"
    #step=3600 # 60*60
    #time_array= np.arange(10800,216000+step,step)
    #npyfilespath ="/home/anthonys/Documents/rico_area_gmean"
    ### becarful here
    #np.save(npyfilespath+str(time_array[file1_numb])+'.npy',area_z_ratio)
    
    #np.save(npyfilespath+str(time_array[20])+'.npy',area_z_ratio)
    ### need to fix time_array and saving

###############################################
end= time.time()
print('Run Time in Seconds:', end-start)










