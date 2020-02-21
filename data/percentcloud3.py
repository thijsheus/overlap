#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 11:41:26 2019

@author: anthonys
"""

"""
calculates overlap ratio base on ht of cld from chopping the cld at the bottom


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

####################################
    
# script

bomexfilenames=[bomextrack36]
ricofilenames=[ricotrack2016]
armfilenames=[armtrack522]
conditional_height=1500


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
            #conditional_height=1500
            index_anvil=np.where(ht > conditional_height)
            index_anvil_size=index_anvil[0].size
            ht_anvil=ht[index_anvil[0]]
            overlap_ratio_anvil=overlap_ratio[index_anvil[0]];
            area_proj_anvil=area_proj[index_anvil[0]]
            print('Clouds that satisfy the given condition: ',index_anvil_size)
            
           
            """
            cbnew1=[];ctnew1=[];cbnew2=[];ctnew2=[];area_projnew1=[];area_projnew2=[];
            overlap_ratio_new1=[];overlap_ratio_new2=[];htnew1=[];htnew2=[];
            index1=[];index2=[];
            
            for i in range(ht.size):
                if ht[i] > conditional_height: #ht[i]*overlap_ratio[i] > 200:
                    cbnew1.append(cb[i])
                    ctnew1.append(ct[i])
                    area_projnew1.append(area_proj[i])
                    overlap_ratio_new1.append(overlap_ratio[i])
                    htnew1.append(ht[i])
                    index1.append(i)
                else:
                    cbnew2.append(cb[i])
                    ctnew2.append(ct[i])
                    area_projnew2.append(area_proj[i])
                    overlap_ratio_new2.append(overlap_ratio[i])
                    htnew2.append(ht[i])
                    index2.append(i)
            
            cbnew1=np.asarray(cbnew1)        
            ctnew1=np.asarray(ctnew1)
            cbnew2=np.asarray(cbnew2)
            ctnew2=np.asarray(ctnew2)
            area_projnew1=np.asarray(area_projnew1)
            area_projnew2=np.asarray(area_projnew2)
            overlap_ratio_new1=np.asarray(overlap_ratio_new1)
            overlap_ratio_new2=np.asarray(overlap_ratio_new2)
            htnew1=np.asarray(htnew1)
            htnew2=np.asarray(htnew2)
            """
            

            
            #overlap_cal, volume_cal, height_cal, projected_area_cal, height_diff = overlap_calculation.overlap_calculation(index1,nrcloudarray,xt,yt,zt)

            
            
            #plt.figure()
            #plt.plot(overlap_cal,height_cal,'o')
            #plt.figure()
            #plt.plot(overlap_ratio_new1,htnew1,'o')



############################################################################
            # calculate overlap of a cld that at some of its height chopped off at the bottom
            
            
            index7=index_anvil[0]
            dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
            gridarea=dx*dy
            gridvol=dx*dy*dz
            nx=xt.size;ny=yt.size;nz=zt.size;
            
            # volume
            #unique_elements, counts_elements = np.unique(nrcloudarray, return_counts=True) # takes cloud # and how many cells, comment out saved 2 sec
            #uni_elmts=unique_elements[1:];count_elmts=counts_elements[1:]  # exclude zero
            # dont need above as we dont need the for loops below
            """
            uninew=[];countnew=[]
            for uni in range(uni_elmts.size):
                for i in range(index7.size):
                    if uni_elmts[uni] == index7[i]: # taking index of cloud w/ certain condition
                        uninew.append(index7[i])
                        countnew.append(count_elmts[index7[i]])
            uninew1=np.asarray(uninew)
            countnew1=np.asarray(countnew)
            """
            cloud_numb=index7 +1 # index is off by 1 as array starts with zero
            #countnew1=count_elmts[index7] # dont need as vol is calculatd in for loop
            
            
            # percent of cloud ht should change countnew1=> change vol, ht, proj area
            #volume_c=dx*dy*dz*countnew1  # calculated later
            
            #cloud_numb=np.array([533]) # choose a particular cld
            #vert_res=25
            height_c=np.zeros(cloud_numb.size)
            overlap_c=np.zeros(cloud_numb.size)
            projected_area_c=np.zeros(cloud_numb.size)
            volume_c=np.zeros(cloud_numb.size)
            
            m=0;
            #alpha=[1.0] #, 0.9, 0.8, 0.7, 0.6, 0.5,0.2]
            #alpha=np.asarray(alpha)
            alpha=np.arange(1,0.15,-0.1) # fraction
            
            #Dcurve=np.zeros((cloud_numb.size,alpha.size))  # multi dimensional array
            for e1 in cloud_numb:  # may take a while to run 
                
                
                location=np.where(nrcloudarray == e1)     # location of all cells with cloud # e1
                locz_size=location[0].size
                #print('locz',locz_size)
                height_c[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
                #projected_area_c[m]=dx*dy*len(np.unique(location[2]+location[1]*nx)) #len of unique # of cells that get proj
                partial_heights=alpha*height_c[m] # fraction of cld height
                #height_diff=height_c[m]-partial_heights
                
                overlap_changed=np.zeros(partial_heights.size)
                projected_area_changed=np.zeros(partial_heights.size)
                vol_changed=np.zeros(partial_heights.size)
                part_loc= -1
                for partial in partial_heights:
                    
                    part_loc= part_loc +1
                    for i in range(locz_size-1,-1,-1):
                        if dz*(-location[0][i]+location[0][-1]+1) <= partial: # finding all cells above chopped cld ht 
                            celltop=i  # celltop is more like cellbottom
                        else:
                            break
                    #if partial - dz*(location[0][celltop]-location[0][0]+1) >0: 
                    difference = partial - dz*(-location[0][celltop]+location[0][-1]+1) # diff btw chopped cld ht and ht of clds w/ full cells
                    #print(difference == 0 )
                    #print(difference)
                    for j in range(celltop,-1,-1):
                        if  location[0][celltop] -1  == location[0][j]: # to use partial cells
                            cellabove=j # cellabove is more like cellbelow
                        elif location[0][celltop] -1  > location[0][j]:
                            break
                    #print('top',celltop)
                    #print('above',cellabove)
                    #loc_loc0=-1
                    
                    #####################################
                    """
                    proj_area_loop1=[]
                    vol_loop1=[]
                    for loc0 in range(locz_size-1,-1,-1):
                        
                        #loc_loc0 = loc_loc0 +1
                        if loc0  >= celltop : 
                            proj_area_loop1.append(location[2][loc0]+location[1][loc0]*nx)
                        elif celltop > loc0 and loc0 >=cellabove:
                            proj_area_loop1.append(location[2][loc0]+location[1][loc0]*nx)
                        if loc0  >= celltop: # full cell vol
                            vol_loop1.append(dx*dy*dz)   
                        elif celltop > loc0 and loc0 >=cellabove: # partial cell vol
                            vol_loop1.append(dx*dy*difference) 
                        else:
                            break
                        
                    proj_area_loop1=np.asarray(proj_area_loop1)
                    #print('loop1',proj_area_loop1.size)
                    vol_loop1=np.asarray(vol_loop1) 
                    
                    ##########################
                    """
                    if partial == height_c[m]:
                        proj_area_loop1=np.zeros(locz_size - celltop)
                        vol_loop1=np.zeros(locz_size  - celltop)
                    else:
                        proj_area_loop1=np.zeros(locz_size -cellabove )
                        vol_loop1=np.zeros(locz_size -cellabove  )
                    
                    loc0_m= -1
                    for loc0 in range(locz_size-1,-1,-1):
                        loc0_m= loc0_m+1
                        #loc_loc0 = loc_loc0 +1
                        if loc0  >= celltop : 
                            proj_area_loop1[loc0_m]=(location[2][loc0]+location[1][loc0]*nx)
                        elif celltop > loc0 and loc0 >=cellabove:
                            proj_area_loop1[loc0_m]=(location[2][loc0]+location[1][loc0]*nx)
                        if loc0  >= celltop: # full cell vol
                            vol_loop1[loc0_m]=(dx*dy*dz)   
                        elif celltop > loc0 and loc0 >=cellabove: # partial cell vol
                            vol_loop1[loc0_m]=(dx*dy*difference) 
                        else:
                            break
                       
                    
                    
                    ######################
                        #if dz*loc0 < parital:
                    #proj_area_loop1=np.asarray(proj_area_loop1)
                    #print('loop1',proj_area_loop1.size)
                    projected_area_changed[part_loc]=dx*dy*len(np.unique(proj_area_loop1))
                    #vol_loop1=np.asarray(vol_loop1)  
                    vol_changed[part_loc]=sum(vol_loop1)
                    overlap_changed[part_loc]=vol_changed[part_loc]/(partial_heights[part_loc]*projected_area_changed[part_loc])
                    if partial == height_c[m]:
                        overlap_c[m]=overlap_changed[part_loc] # make sure calculated overlap is correct for full ht
                        volume_c[m]=vol_changed[part_loc]
                        projected_area_c[m]=projected_area_changed[part_loc]
                        
                #plt.figure()
                #if e1==351:
                
                overlap_changed_grad=np.gradient(overlap_changed)
                alpha_grad=np.gradient(alpha)
                
                #Dcurve[m]=abs(alpha_grad/overlap_changed_grad)
                
                plt.plot(overlap_changed,partial_heights/height_c[m],'o-',label=str(cloud_numb[m]))
                plt.title('fractional height (cut at bottom) vs. overlap')
                
                plt.legend(loc='upper right', bbox_to_anchor=(1.40, 1.0), fontsize='small',ncol=2)
                
                
                
                
                m = m+1
            
            
            
            overlap_c=volume_c/(height_c*projected_area_c)


















###################################

end= time.time()
print('Run Time in Seconds:', end-start)



















