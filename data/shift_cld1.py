#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 12:35:46 2019

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
import cProfile
from numpy import percentile
#import overlap_calculation
import sys

from scipy.stats import spearmanr


start=time.time()



###############################################################





epsilon1=50     ##### helps find clds that are split by grid box
#arbitrarily set, most outliers are > 200, most nonoutliers clds are about 0

#file1_numb=-1

"""
datau = Dataset('/data/rico/u.nc','r')
#u=datau.variables['u'][:,:,:,:]
datav = Dataset('/data/rico/v.nc','r')
#v=datav.variables['v'][:,:,:,:]
dataw = Dataset('/data/rico/w.nc','r')
time_t = datau.variables['time'][:]
step=3600;time_array= np.arange(step,216000+step,step)
case='rico'
"""
"""
datau = Dataset('/data/bomex/u.nc','r')
#u=datau.variables['u'][:,:,:,:]
datav = Dataset('/data/bomex/v.nc','r')
#v=datav.variables['v'][:,:,:,:]
dataw = Dataset('/data/bomex/w.nc','r')
time_t = datau.variables['time'][:]
step=1800;time_array= np.arange(step,36000+step,step)
case='bomex'
"""
"""
datau = Dataset('/data/lasso/sims/20160830/u.nc','r')
#u=datau.variables['u'][:,:,:,:]
datav = Dataset('/data/lasso/sims/20160830/v.nc','r')
#v=datav.variables['v'][:,:,:,:]
dataw = Dataset('/data/lasso/sims/20160830/w.nc','r')
"""

datau = Dataset('/data/arm/u.nc','r')
#u=datau.variables['u'][:,:,:,:]
datav = Dataset('/data/arm/v.nc','r')
#v=datav.variables['v'][:,:,:,:]
dataw = Dataset('/data/arm/w.nc','r')
time_t = datau.variables['time'][:]
step=1800;time_array= np.arange(10800,52200+step,step)
case='arm'


conditional_height=50
#file1_numb=-1
begin=5
### accessing multiple datasets easily
#Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
#Bfilenames = Afilenames[begin:]
#file3_numb = np.arange(5,19+1)
#Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
#begin=22
#Bfilenames = Afilenames[22:23]
#Bfilenames = Afilenames[55:56]
#Bfilenames = [Afilenames[22], Afilenames[55]]
#Bfilenames = [Afilenames[22]]
#file1_numb=begin-1 ## python: arrays start at index 0
#file3_numb = [22,55]
#file3_numb = [22]
Afilenames = sorted(glob.glob('/data/arm/*.track.nc'))
#Bfilenames = Afilenames[1:2]
Bfilenames = [Afilenames[10], Afilenames[17]]
file3_numb = [10,17]

file2_numb = -1


#############################
#### loop over trackfiles
for file1 in Bfilenames:
    print(file1)

    data = Dataset(file1,'r')

    file2_numb = file2_numb+1
    file1_numb = file3_numb[file2_numb]
    
    
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
    
    #index1=np.where( (400 < ht) & (ht < 500)) 
    index1=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
    index1_size=index1[0].size
    ht=ht[index1[0]]   # taking the ht values according to indices above
    overlap_ratio=overlap_ratio[index1[0]];
    area_proj=area_proj[index1[0]]
    print('Clouds that satisfy the given condition: ',index1_size)
    print('conditional height is',conditional_height,'meters')


    cloud_numb=index1[0] +1 # index is off by 1 as array starts with zero
            
    height_c=np.zeros(cloud_numb.size)
    overlap_c=np.zeros(cloud_numb.size)
    projected_area_c=np.zeros(cloud_numb.size)
    volume_c=np.zeros(cloud_numb.size)
    m=-1;
    
    alpha=np.arange(1,0,-1) # just 1
    
    
    
    overlap_changed_matrix=np.zeros((cloud_numb.size,alpha.size))
    area_z_ratio=np.zeros(cloud_numb.size)
    shift_max=np.zeros(cloud_numb.size)
    shift_min=np.zeros(cloud_numb.size)
    shift_avg=np.zeros(cloud_numb.size)
    shift_distance=np.zeros(cloud_numb.size)
    
    for e1 in cloud_numb:  # may take a while to run 
        
        m=m+1
        location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1
        layers=(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z 
        #layers=height/dz
        ###
        
        base=np.amin(location[0]);top=np.amax(location[0]);
        #levels=range(base,top+1,1)
        
        x_unique=np.unique(location[2])
        y_unique=np.unique(location[1])
        xtest1=abs(np.mean(x_unique)-np.median(x_unique))  # median is resistant while mean is not
        ytest1=abs(np.mean(y_unique)-np.median(y_unique))
 
        #xtest1=abs(np.mean(location[2])-np.median(location[2]))
        #ytest1=abs(np.mean(location[1])-np.median(location[1]))
        #x_small=np.extract(location[2] < 512,location[2])
        #x_large=np.extract(x_unique > 512,x_unique)
        #y_small=np.extract(y_unique < 512,y_unique)
        #y_large=np.extract(y_unique > 512,y_unique)
        
        if xtest1 > epsilon1:
            splitx=np.where(location[2] < 512)
            location[2][splitx] = location[2][splitx] + nx
            
        if ytest1 > epsilon1:
            splity=np.where(location[1] < 512)
            location[1][splity] = location[1][splity] + ny
        
        
        #sys.exit("Error message")
        

        location_matrix=np.zeros((location[0].size,3))
        location_matrix[:,0]=location[0] # z coor
        location_matrix[:,1]=location[1] # y coor
        location_matrix[:,2]=location[2] # x coor
        COM=np.zeros((top-base+1,3))
        newlocation_matrix=np.zeros(location_matrix.shape)
        cross_area_z=np.zeros(int(layers))
        #u, indices = np.unique(location[0], return_index=True)
        #COM[:,0]=u
        
        
        for z in np.arange(base,top+1,1):
            findz=np.where(location[0] == z)
            cross_area_z[z-base]=dx*dy*len(findz[0])
            # height / dz = layers
            
            
            #COM[z-base,:]=[z,sum(location[1][findz[0]])/findz[0].size, sum(location[2][findz[0]])/findz[0].size]
            COM[z-base,:]=[z,np.mean(location[1][findz[0]]), np.mean(location[2][findz[0]])] #find center of mass w/ arithmetic mean
            ### not quite (close)
            """
            if z<top:
                COM[z-base,1] = sum(location[1][indices[z-base]: indices[z-base+1]-1]) / len(location[1][indices[z-base]: indices[z-base+1]-1])
                COM[z-base,2] = sum(location[2][indices[z-base]: indices[z-base+1]-1]) / len(location[2][indices[z-base]: indices[z-base+1]-1])
            else:
                tempz=np.where(location[0] == z)
                COM[z-base,:]=[z,sum(location[1][tempz[0]])/len(tempz[0]), sum(location[2][tempz[0]])/len(tempz[0])]
        
            """
        area_z_ratio[m]= np.mean(cross_area_z)/np.amax(cross_area_z)    
        center_choice = [0,0]
        #center_choice=[np.amin(COM[:,1]), np.amin(COM[:,2])] 
        #center_choice=[np.mean(COM[:,1]), np.mean(COM[:,2])] 
        #movement=np.subtract(center_choice , COM[:,1:])  #### need multiply by dx & dy
        movement = 25*COM[:,1:] ### having center at origin,spacing in x and y direction
        
        
        shift_cld= np.sqrt(movement[:,0]**2 + movement[:,1]**2)
        #shift_cld = (abs(movement[:,0])  + abs(movement[:,1]) ) / 2
        shift_max[m] = np.amax(shift_cld)
        shift_min[m] = np.amin(shift_cld)
        shift_avg[m] = np.mean(shift_cld)
        shift_distance[m] = np.amax(shift_cld) - np.amin(shift_cld) 
        
    """
    #### save data
    npyfilespath ="/home/anthonys/Documents/"
    
    name1='shift_max_'
    name2='shift_min_'
    name3='shift_avg_'
    name4='shift_distance'
    #name5='area_z_ratio_'
    
    
    ### 
    np.save(npyfilespath+ name1  + case +str(time_array[file1_numb])+'.npy',shift_max)
    np.save(npyfilespath+ name2  + case +str(time_array[file1_numb])+'.npy',shift_min)
    np.save(npyfilespath+ name3  + case +str(time_array[file1_numb])+'.npy',shift_avg)
    np.save(npyfilespath+ name4  + case +str(time_array[file1_numb])+'.npy',shift_distance)
    #np.save(npyfilespath+ name5  + case +str(time_array[file1_numb])+'.npy',area_z_ratio)
    """ 




###################################

end= time.time()
print('Run Time in Seconds:', end-start)





























