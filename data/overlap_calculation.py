#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 11:50:04 2019

@author: anthonys
"""

#from netCDF4 import Dataset
import numpy as np 
#import struct 
#import netCDF4
from netCDF4 import Dataset
import netCDF4 as nc
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
import collections
from collections import Counter
##############################################


"""
def overlap_calculation(index,nrcloudarray,xt,yt,zt):

    
    index7=np.asarray(index)
    dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
    gridarea=dx*dy
    gridvol=dx*dy*dz
    nx=xt.size;ny=yt.size;nz=zt.size;
    
    # volume
    unique_elements, counts_elements = np.unique(nrcloudarray, return_counts=True) # takes cloud # and how many cells
    uni_elmts=unique_elements[1:];count_elmts=counts_elements[1:]  # exclude zero
    uninew=[];countnew=[]
    for uni in range(uni_elmts.size):
        for i in range(index7.size):
            if uni_elmts[uni] == index7[i]:
                uninew.append(index7[i])
                countnew.append(count_elmts[index7[i]])
    uninew1=np.asarray(uninew)
    countnew1=np.asarray(countnew)
    
    # percent of cloud ht should change countnew1=> change vol, ht, proj area
    volume_c=dx*dy*dz*countnew1
    
    
    #vert_res=25
    height_c=np.zeros(uninew1.size)
    projected_area_c=np.zeros(uninew1.size)
    m=0;
    for e1 in uninew1:  # may take a while to run 
        location=np.where(nrcloudarray == e1+1)    # location of all cells with cloud # e1+1, e1 is the index and cloud # is index + 1 
        height_c[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
        projected_area_c[m]=dx*dy*len(np.unique(location[2]+location[1]*nx)) #len of unique # of cells that get proj
        m = m+1
    
    
    
    overlap_c=volume_c/(height_c*projected_area_c)

    return overlap_c, volume_c, height_c, projected_area_c
"""

def overlap_calculation(index1,nrcloudarray,xt,yt,zt):

    
    index7=np.asarray(index1)
    dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
    gridarea=dx*dy
    gridvol=dx*dy*dz
    nx=xt.size;ny=yt.size;nz=zt.size;
    
    # volume
    unique_elements, counts_elements = np.unique(nrcloudarray, return_counts=True) # takes cloud # and how many cells
    uni_elmts=unique_elements[1:];count_elmts=counts_elements[1:]  # exclude zero
    uninew=[];countnew=[]
    for uni in range(uni_elmts.size):
        for i in range(index7.size):
            if uni_elmts[uni] == index7[i]:
                uninew.append(index7[i])
                countnew.append(count_elmts[index7[i]])
    uninew1=np.asarray(uninew)
    countnew1=np.asarray(countnew)
    
    # percent of cloud ht should change countnew1=> change vol, ht, proj area
    volume_c=dx*dy*dz*countnew1
    
    
    #vert_res=25
    height_c=np.zeros(uninew1.size)
    projected_area_c=np.zeros(uninew1.size)
    m=0;
    #alpha=[1.0] #, 0.9, 0.8, 0.7, 0.6, 0.5,0.2]
    #alpha=np.asarray(alpha)
    alpha=np.arange(1,0.15,-0.05)
    for e1 in uninew1:  # may take a while to run 
        
        
        location=np.where(nrcloudarray == e1+1)    # location of all cells with cloud # e1+1, e1 is the index and cloud # is index + 1 
        height_c[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
        projected_area_c[m]=dx*dy*len(np.unique(location[2]+location[1]*nx)) #len of unique # of cells that get proj
        partial_heights=alpha*height_c[m]
        height_diff=height_c[m]-partial_heights
        
        overlap_changed=np.zeros(partial_heights.size)
        projected_area_changed=np.zeros(partial_heights.size)
        vol_changed=np.zeros(partial_heights.size)
        part_loc= -1
        for partial in partial_heights:
            celltop=0
            part_loc= part_loc +1
            proj_area_loop1=[]
            vol_loop1=[]
            for i in range(location[0].size):
                if dz*(location[0][i]-location[0][0]+1) <= partial:
                    celltop=i
            if partial - dz*(location[0][celltop]-location[0][0]+1) >0:
                difference = partial - dz*(location[0][celltop]-location[0][0]+1)
                for j in range(celltop,location[0].size):
                    if  location[0][celltop] +1  == location[0][j]:
                        cellabove=j
            #loc_loc0=-1
            for loc0 in range(location[0].size):
                
                #loc_loc0 = loc_loc0 +1
                if loc0  <= celltop :
                    proj_area_loop1.append(location[2][loc0]+location[1][loc0]*nx)
                if celltop < loc0 and loc0 <=cellabove:
                    proj_area_loop1.append(location[2][loc0]+location[1][loc0]*nx)
                if loc0  <= celltop:
                    vol_loop1.append(dx*dy*dz)   
                if celltop < loc0 and loc0 <=cellabove:
                    vol_loop1.append(dx*dy*difference) 
                #if dz*loc0 < parital:
            proj_area_loop1=np.asarray(proj_area_loop1)
            projected_area_changed[part_loc]=dx*dy*len(np.unique(proj_area_loop1))
            vol_loop1=np.asarray(vol_loop1)  
            vol_changed[part_loc]=sum(vol_loop1)
            overlap_changed[part_loc]=vol_changed[part_loc]/(partial_heights[part_loc]*projected_area_changed[part_loc])
        plt.plot(overlap_changed,alpha,'o-')
        #plt.xlim([0, 1])
        #plt.plot(overlap_changed,height_c[m],'o')
        #plt.plot(overlap_ratio_new1,htnew1,'o')
        m = m+1
    
    
    
    overlap_c=volume_c/(height_c*projected_area_c)


    
    
    

    return overlap_c, volume_c, height_c, projected_area_c, height_diff

































