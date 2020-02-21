#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:18:45 2019

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

#import overlap_calculation

start=time.time()


conditional_height=100

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
    
    # selecting cld to plot explicitly
    #cloud_numb=np.array([278,313,351]) #ricotrack2016
    #cloud_numb=np.array([137,371,431]) #ricotrack2016
    #cloud_numb=np.array([70,254,486])   #ricotrack2124
    #cloud_numb=np.array([196,199,242])   #ricotrack2124
    #cloud_numb=np.array([533]) # ricotrack828
    #cloud_numb=np.array([4320, 4526]) # ricotrack2016
    #cloud_numb=np.array([224]) # ricotrack900
    #cloud_numb=np.array([330]) # ricotrack1548
    #cloud_numb=np.array([575]) # ricotrack828
    
    
    height_convex=np.zeros(cloud_numb.size)
    vol_convex=np.zeros(cloud_numb.size)
    projarea_convex=np.zeros(cloud_numb.size)
    overlap_convex=np.zeros(cloud_numb.size)
    ######################################################## convex hull
    
    m=0;
    for e1 in cloud_numb:  # may take a while to run 
        
        location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1
        location_matrix=np.zeros((location[0].size,3))
        location_matrix[:,0]=location[0]
        location_matrix[:,1]=location[1]
        location_matrix[:,2]=location[2]
        hull_3d=ConvexHull(location_matrix, qhull_options='QJ')
        height_convex[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
        vol_convex[m]=dx*dy*dz*hull_3d.volume
        bd_pts=hull_3d.points[hull_3d.vertices]
        hull_2d=ConvexHull(bd_pts[:,1:], qhull_options='QJ')     
        projarea_convex[m]=dx*dy*hull_2d.volume   # volume and area are named base in 3d, so the 2d volume is indeed area
        overlap_convex[m]= vol_convex[m] / (height_convex[m]*projarea_convex[m])
        m=m+1
    
    overlap_contra=1/overlap_ratio - 1/overlap_convex  
    
    step=1800 # 30*60
    time_array= np.arange(10800,36000+step,step)
    npyfilespath ="/home/anthonys/Documents/bomex_overlap_convex"
    #np.save(npyfilespath+str(time_array[file1_numb])+'.npy',overlap_convex)
    
    ### plotting in convex hull in 3d and 2d
    """
    fig=plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(dx*bd_pts[:,2], dy*bd_pts[:,1],dz*bd_pts[:,0], color='green');
    #ax.scatter3D(dx*location_matrix[:,2], dy*location_matrix[:,1],dz*location_matrix[:,0], color='blue');
    #ax.plot_trisurf(dx*bd_pts[:,2], dy*bd_pts[:,1],dz*bd_pts[:,0], color='grey'); 
    for simplex in hull_3d.simplices:
        #fig=plt.figure()
        #ax = plt.axes(projection='3d')
        
        ax.plot3D(dx*location_matrix[simplex,2], dy*location_matrix[simplex,1],dz*location_matrix[simplex,0], color='black'); 
        #ax.plot_trisurf(dx*location_matrix[simplex,2], dy*location_matrix[simplex,1],dz*location_matrix[simplex,0], color='grey'); 
        #plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
    
    plt.title('3d Convex Hull')
    plt.figure()   
    plt.plot(dx*bd_pts[:,2], dy*bd_pts[:,1], 'o')
    for simplex in hull_2d.simplices:
        plt.plot(dx*bd_pts[simplex, 2], dy*bd_pts[simplex, 1], 'k-')
    plt.title('2d Convex Hull')
    """
    
    """
    ##########################################################
    ### plot overlap
    bins=dz
    
    plt.figure()
    plt.hist2d(overlap_convex,height_convex,bins=bins,cmin=0.5)
    plt.title('Height vs. Convex Hull Overlap ')
    plt.xlim([0,1])
    colorbar = plt.colorbar()
    colorbar.set_label('counts in bin')
    colorbar.set_label('counts in bin')
    
    
    plt.figure()
    plt.hist2d(overlap_ratio,ht,bins=bins,cmin=0.5)
    plt.title('Height vs. Actual Overlap')
    plt.xlim([0,1])
    colorbar = plt.colorbar()
    colorbar.set_label('counts in bin')
    """
    
    
end= time.time()
print('Run Time in Seconds:', end-start)




















