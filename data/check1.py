#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 10:36:00 2019

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
import os
import sys
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import time
#import pkgutil
#import collections
#from collections import Counter
from scipy.spatial import ConvexHull
import cProfile
from numpy import percentile
#import overlap_calculation
from scipy.stats import spearmanr


start=time.time()


#############################################################################

# definitions
def pyth(u,v):    # magnitude
    return np.sqrt(u*u+v*v)
    



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
ricotrack144 = Dataset('/data/rico/l.0014400.track.nc','r')
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


filenames=[bomexd, ricod, armd]

bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72, bomextrack90, bomextrack108, bomextrack126, bomextrack144, bomextrack162, bomextrack180, bomextrack198, bomextrack216]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144, ricotrack180, ricotrack216, ricotrack252, ricotrack288, ricotrack324, ricotrack360, ricotrack396]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162, armtrack180, armtrack198, armtrack216, armtrack234, armtrack252, armtrack270, armtrack288]

###########################################################################

# script
filenames=[ricod]
ricofilenames=[ricotrack828]
bomexfilenames=[bomextrack360]
armfilenames=[armtrack522]
#conditional_height=1000
epsilon1=50    #arbitrarily set, most outliers are > 200, most nonoutliers clds are about 0

kt=0;
file_numb=-1
file1_numb=-1

for file in filenames:
    file_numb = file_numb+1
    #zt=file.variables['z'][:]
    zh=file.variables['zh'][:]
    time_t=file.variables['time'][:]
    u=file.variables['u'][:,:]
    v=file.variables['v'][:,:]
    w=file.variables['w'][:,:]
    
    
            
    if file == ricod:
        for file1 in ricofilenames:
            file1_numb = file1_numb+1
            ht=file1.variables['ht'][:]
            cb=file1.variables['cb'][:]
            ct=file1.variables['ct'][:]
            cv=file1.variables['cv'][:]
            cp=file1.variables['cp'][:]
            overlap_ratio=file1.variables['chr'][:]
            area_proj=file1.variables['area_proj'][:]
            nrcloud=file1.variables['nrcloud'][:,:,:]
            cfrac=file1.variables['cfrac'][:]
            zt=file1.variables['z'][:]
            xt=file1.variables['x'][:]
            yt=file1.variables['y'][:]
            nr=file1.variables['nr'][:]
            cld_mask=file1.variables['cld_mask'][:,:,:]
            
            
            nrcloudarray = np.ma.getdata(nrcloud) # unmask array
            
            
            
            condition_vec=[100]#, 500, 100, 50] # to go through multiple conditions quickly
            
            for conditional_height in condition_vec:
                
                #conditional_height=1500
                index_greater=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
                index_greater_size=index_greater[0].size
                ht_anvil=ht[index_greater[0]]   # taking the ht values according to indices above
                overlap_ratio_anvil=overlap_ratio[index_greater[0]];
                area_proj_anvil=area_proj[index_greater[0]]
                
                print('Total Number of Clouds in dataset is %.0f' % nr.size)
                print('Condition: Clouds with height greater than',conditional_height,'meters')
                print('Clouds that satisfy the given condition: ',index_greater_size)




"""
### combines npy files into one
### only use once, affect will be compounded, delete combined npy file then run again 
fpath ="path_Of_my_final_Big_File"
npyfilespath ="/home/anthonys/Documents/Temp1"   
os.chdir(npyfilespath)
npfiles= glob.glob("*.npy")
npfiles.sort()
all_arrays = []
for i, npfile in enumerate(npfiles):
    all_arrays.append(np.load(os.path.join(npyfilespath, npfile)))
np.save(fpath, np.concatenate(all_arrays,axis=0))
"""
######

"""
### accessing multiple datasets easily
Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
Bfilenames = Afilenames[5:]
for f in Bfilenames:
    print(f)

    data = Dataset(f,'r')
step=1800 # 30*60
time_array= np.arange(10800,36000+step,step)
"""
"""
######################################################
### Concatenate data
### accessing multiple datasets easily
Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
Bfilenames = Afilenames[5:]
arrays1=[];arrays2=[];arrays3=[]
for f in Bfilenames:
    print(f)

    data = Dataset(f,'r')
    #overlap_ratio=data.variables['chr'][:]
    #arrays1.append(overlap_ratio)
    #ht=data.variables['ht'][:]
    #arrays2.append(ht)
    cp=data.variables['cp'][:]
    arrays3.append(cp)

#Bomex_overlap_ratio=np.concatenate(arrays1,axis=0)
#Bomex_ht=np.concatenate(arrays2,axis=0)
Bomex_cp=np.concatenate(arrays3,axis=0)
"""
########
"""
#### accessing multiple npy files easily
Afilenames = sorted(glob.glob('/home/anthonys/Documents/overlap_no_shear_bomex*.npy'))
#Bfilenames = Afilenames[5:]
Bfilenames = [x for i,x in enumerate(Afilenames) if i!=14] # exclude 14th elt
arrays=[]
for f in Bfilenames:
    print(f)

    data = np.load(f)
    arrays.append(data)
    #print(data.size)
Bomex_shear=np.concatenate(arrays,axis=0)

#######

#### accessing multiple npy files easily
Afilenames = sorted(glob.glob('/home/anthonys/Documents/area_z_ratio_bomex*.npy'))
#Bfilenames = Afilenames[5:]
Bfilenames = [x for i,x in enumerate(Afilenames) if i!=14] # exclude 14th elt
arrays=[]
for f in Bfilenames:
    print(f)

    data = np.load(f)
    arrays.append(data)
    #print(data.size)
Bomex_areaz=np.concatenate(arrays,axis=0)
"""
"""
#### accessing multiple npy files easily
Afilenames = sorted(glob.glob('/home/anthonys/Documents/bomex_overlap_convex*.npy'))
#Bfilenames = Afilenames[5:]
#Bfilenames = [x for i,x in enumerate(Afilenames) if i!=14] # exclude 14th elt
Bfilenames = Afilenames
arrays=[]
for f in Bfilenames:
    print(f)

    data = np.load(f)
    arrays.append(data)
    #print(data.size)
Bomex_overlap_convex=np.concatenate(arrays,axis=0)
"""

"""
#### accessing multiple npy files easily
Afilenames = sorted(glob.glob('/home/anthonys/Documents/shift_avg_bomex*.npy'))
Bfilenames = Afilenames[5:]
#Bfilenames = [x for i,x in enumerate(Afilenames) if i!=14] # exclude 14th elt
Bfilenames = Afilenames
arrays=[]
for f in Bfilenames:
    print(f)

    data = np.load(f)
    arrays.append(data)
    #print(data.size)
Bomex_shift_avg=np.concatenate(arrays,axis=0)
"""
#### need to save concatenated data
#######
########################################################################


##############################
####################################################
###############################################
end= time.time()
print('Run Time in Seconds:', end-start)
