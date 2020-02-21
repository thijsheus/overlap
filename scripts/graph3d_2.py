#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:13:44 2019

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
#import glob
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
# def 

def overlap(s,h,l):
    #return l / (l +s*h)
    return 1 - (s*h / (l + s*h))
####################################
    
# script
filenames=[ricod]
bomexfilenames=[bomextrack36]
ricofilenames=[ricotrack828]
armfilenames=[armtrack522]
conditional_height=2000

### choose which case and track file


for file in filenames:
    #zt=file.variables['z'][:]
    zh=file.variables['zh'][:]
    time_t=file.variables['time'][:]
    u=file.variables['u'][:,:]
    v=file.variables['v'][:,:]
    w=file.variables['w'][:,:]

    
            
    if file == ricod:
        for file1 in ricofilenames:
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
            
            nrcloudarray = np.ma.getdata(nrcloud)
            #conditional_height=2000
            index_anvil=np.where(ht > conditional_height)
            index_anvil_size=index_anvil[0].size
            ht_anvil=ht[index_anvil[0]]
            overlap_ratio_anvil=overlap_ratio[index_anvil[0]];
            area_proj_anvil=area_proj[index_anvil[0]]
            print('Clouds that satisfy the given condition: ',index_anvil_size)
            


            index7=index_anvil[0]
            dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
            gridarea=dx*dy
            gridvol=dx*dy*dz
            nx=xt.size;ny=yt.size;nz=zt.size;
            
            
            cloud_numb=index7 +1 # index is off by 1 as array starts with zero
            
            # selecting cld to plot explicitly
            #cloud_numb=np.array([278,313,351]) #ricotrack2016
            #cloud_numb=np.array([137,371,431]) #ricotrack2016
            #cloud_numb=np.array([70,254,486])   #ricotrack2124
            #cloud_numb=np.array([196,199,242])   #ricotrack2124
            #cloud_numb=np.array([533]) # ricotrack828
            #cloud_numb=np.array([4320, 4526]) # ricotrack2016
            #cloud_numb=np.array([224]) # ricotrack900
            #cloud_numb=np.array([330]) # ricotrack1548
            cloud_numb=np.array([533,860,1991,2004]) # ricotrack828 
            
           
          
            m=0
            for e1 in cloud_numb:  # may take a while to run 
                location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1+1, e1 is the index and cloud # is index + 1 
                
                
                
                """
                ### plotting in cld in 3d
                
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.scatter3D(dx*location[2], dy*location[1], dz*location[0], c=dz*location[0] ,cmap='Paired'); 
                #cmap= single color:copper, cool, winter, multicolor:  Dark2, Paired
                #ax.plot_trisurf(location[2], location[1], location[0], cmap='viridis', edgecolor='none');
                plt.title(str(cloud_numb[m]))
                plt.show()
                """
                
                #### ploting z vs. x or z vs. x in scatter plots
                """
                plt.figure()
                plt.plot(dx*location[2], dz*location[0],'go') 
                plt.xlabel('x');plt.ylabel('z')
                plt.title(str(cloud_numb[m]))
                plt.figure()
                plt.plot(dy*location[1], dz*location[0],'go')
                plt.xlabel('y');plt.ylabel('z')
                plt.title(str(cloud_numb[m]))
                """
                #### ploting z vs. x or z vs. x in 2d historgrams
                
                bins=dz
                plt.figure()
                plt.hist2d(dx*location[2],dz*location[0],bins=bins,cmin=0.5)
                plt.xlabel('x');plt.ylabel('z')
                plt.title(str(cloud_numb[m]))
                colorbar = plt.colorbar()
                colorbar.set_label('counts in bin')
                plt.figure()
                plt.hist2d(dy*location[1],dz*location[0],bins=bins,cmin=0.5)
                plt.xlabel('y');plt.ylabel('z')
                plt.title(str(cloud_numb[m]))
                colorbar = plt.colorbar()
                colorbar.set_label('counts in bin')
                
                
                m=m+1
            
            
            
            
end= time.time()
print('Run Time in Seconds:', end-start)


