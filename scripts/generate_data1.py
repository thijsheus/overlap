#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:46:53 2019

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


#############################################################################

# definitions
def pyth(u,v):    # magnitude
    return np.sqrt(u*u+v*v)
    



###################################################
conditional_height=0
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

#begin=5
### accessing multiple datasets easily
#Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
#Bfilenames = Afilenames[begin:]
#Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:8]
#Bfilenames = Afilenames[22:23]
#Bfilenames = Afilenames[1:]
#Bfilenames = [Afilenames[22], Afilenames[55]]
#Afilenames = sorted(glob.glob('/data/lasso/sims/20160830/*.track.nc'))
#Bfilenames = Afilenames[1:]
#Bfilenames = [Afilenames[23], Afilenames[41]]
Afilenames = sorted(glob.glob('/data/arm/*.track.nc'))
#Bfilenames = Afilenames[1:2]
Bfilenames = [Afilenames[10], Afilenames[17]]

#file1_numb=begin-1 ## python: arrays start at index 0
file2_numb = -1
file3_numb = [10,17]

################################



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

    ######################################################################
    #### wind and udraft calculations
    
    ###
    access = np.where(time_t == time_array[file1_numb])[0][0]
    u=datau.variables['u'][access,:,:,:]
    xh=datau.variables['xh'][:]
    v=datav.variables['v'][access,:,:,:]
    yh=datav.variables['yh'][:]
    w=dataw.variables['w'][access,:,:,:]
    zh=dataw.variables['zh'][:]
    ###
    
    uint=interp1d(xh,u,axis=2,fill_value='extrapolate')
    vint=interp1d(yh,v,axis=1,fill_value='extrapolate')
    wint=interp1d(zh,w,axis=0,fill_value='extrapolate')
    
    ut=uint(xt)
    vt=vint(yt)
    wt=wint(zt)
    
    #w=interp1d(zh,w,axis=0,fill_value='extrapolate')
    #ut=u;vt=v;#wt=w;
    
    
    ut1=np.sum(np.sum(ut,axis=1),axis=1)/(nx*ny)
    vt1=np.sum(np.sum(vt,axis=1),axis=1)/(nx*ny)
    #ut1=np.mean(np.mean(ut,axis=1),axis=1)
    #vt1=np.mean(np.mean(vt,axis=1),axis=1)
    #wt1=np.mean(np.mean(wt,axis=1),axis=1)
    
    
    uz=ut1
    vz=vt1
    #wz=wt1
    
    #uz=np.mean(u,axis=0) # avg among diff times
    #vz=np.mean(v,axis=0) # avg among diff times
    shear_calu=np.zeros(cb.size)
    shear_calv=np.zeros(cb.size)
    """
    for i in range(cb.size):
        shear_calu[i] = (uz[ct[i]] - uz[cb[i]])
        shear_calv[i]= (vz[ct[i]] - vz[cb[i]])
    """
    
    
    #shear_cal = np.sqrt( shear_calu**2 + shear_calv**2  )
    
    #shear_cal = abs(shear_calu + shear_calv) / 2
    
    
    
    #######
    
    duz=np.gradient(uz) # central diff in uz
    dvz=np.gradient(vz)
    
    
    
    #s_u=abs(duz);s_v=abs(dvz)
    s_u=duz;s_v=dvz;
    
    #Gmean=(s_u*s_v)**(0.5)
    Amean=(s_u+s_v)/2;#Amean=abs(Amean)
    #H=Gmean**2/Amean
    Hmean=(2*s_u*s_v) / (s_u + s_v);#Hmean=abs(Hmean)
    #pyth= (s_u**2 + s_v**2)**0.5
    
    shear0=np.zeros(cb.size)
    shear00=np.zeros(cb.size)
    for i in range(cb.size): ##### calculating shear
        
        shear0[i] = sum(Hmean[cb[i]:ct[i]+1])    
        shear00[i] = sum(Amean[cb[i]:ct[i]+1])
        
        shear_calu[i] = (uz[ct[i]] - uz[cb[i]])
        shear_calv[i]= (vz[ct[i]] - vz[cb[i]])
    
    
    shear_cal = abs(shear_calu + shear_calv) / 2
    
    ###########################################################################

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
    
    ### not chopping so just set alpha=1
    alpha=np.arange(1,0,-1) # just 1

    overlap_changed_matrix=np.zeros((cloud_numb.size,alpha.size))
    area_z_ratio=np.zeros(cloud_numb.size)
    
    
    height_convex=np.zeros(cloud_numb.size)
    vol_convex=np.zeros(cloud_numb.size)
    projarea_convex=np.zeros(cloud_numb.size)
    overlap_convex=np.zeros(cloud_numb.size)
    
    w_avg = np.zeros(index1_size)
    wz_max = np.zeros(index1_size)
    wz_cb = np.zeros(index1_size)
    wz_ct = np.zeros(index1_size)
    wz_avg = np.zeros(index1_size)
    wz_min = np.zeros(index1_size)
    
    shift_max=np.zeros(cloud_numb.size)
    shift_min=np.zeros(cloud_numb.size)
    shift_avg=np.zeros(cloud_numb.size)
    shift_distance=np.zeros(cloud_numb.size)
    
    m=-1; #m=0;
    
    #### loop over clds
    for e1 in cloud_numb:  # may take a while to run 
        
        m = m+1
        
        location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1
        layers=(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z 
        #layers=height/dz
        ###
        
        base=np.amin(location[0]);top=np.amax(location[0]);
        #levels=range(base,top+1,1)
        
        #### test to find clds that are split by the grid box 
        x_unique=np.unique(location[2])
        y_unique=np.unique(location[1])
        xtest1=abs(np.mean(x_unique)-np.median(x_unique))  # median is resistant while mean is not
        ytest1=abs(np.mean(y_unique)-np.median(y_unique))  
 
        #### reconnnects split cld , idea: moves part of cld to otherside
        if xtest1 > epsilon1:
            splitx=np.where(location[2] < nx/2)  
            location[2][splitx] = location[2][splitx] + nx
            
        if ytest1 > epsilon1:
            splity=np.where(location[1] < ny/2)
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
        
        ############################################################
        #### generate updraft speed data
        basez=base;topz=top;basey=min(location[1]);topy=max(location[1]);basex=min(location[2]);topx=max(location[2]);
        w_abs=abs(wt)
        #wz_avg=np.sum(w_abs[basez:topz+1,:,:],axis=0) / (topz - basez+1)
        #wy_avg=np.sum(wz_avg[basey:topy+1,:],axis=0) / (topy - basey+1)
        #wx_avg=np.sum(wy_avg[basex:topx+1],axis=0) / (topx - basex+1)
        ##wz_avg=np.mean(w_abs[min(location[0]):max(location[0])+1,:,:],axis=0)
        ##wy_avg=np.mean(wz_avg[min(location[1]):max(location[1])+1,:],axis=0)
        ##wx_avg=np.mean(wy_avg[min(location[2]):max(location[2])+1],axis=0)
        #w_avg[m]=wx_avg
        #print(w_avg[m])
        """
        wy_avg=np.sum(w_abs[:,basey:topy+1,:],axis=1) / (topy - basey+1)
        wx_avg=np.sum(wy_avg[:,basex:topx+1],axis=1) / (topx - basex+1)
        wz = wx_avg[basez:topz+1]
        
        wz_max[m] = np.amax(wz)
        wz_cb[m] = wz[0]
        wz_ct[m] = wz[-1]
        wz_avg[m] = np.mean(wz)
        wz_min[m] = np.amin(wz)
        """
        #####################
        ### generate Convex hull data
        if dz*layers >100:
            hull_3d=ConvexHull(location_matrix, qhull_options='QJ')
            height_convex[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
            vol_convex[m]=dx*dy*dz*hull_3d.volume
            bd_pts=hull_3d.points[hull_3d.vertices]
            hull_2d=ConvexHull(bd_pts[:,1:], qhull_options='QJ')     
            projarea_convex[m]=dx*dy*hull_2d.volume   # volume and area are named base in 3d, so the 2d volume is indeed area
            overlap_convex[m]= vol_convex[m] / (height_convex[m]*projarea_convex[m])
        
        ########################
        ##### find center of mass and area of each layer of cld 
        for z in np.arange(base,top+1,1):
            findz=np.where(location[0] == z)
            cross_area_z[z-base]=dx*dy*len(findz[0])
            # height / dz = layers
            
            #COM[z-base,:]=[z,sum(location[1][findz[0]])/findz[0].size, sum(location[2][findz[0]])/findz[0].size]
            COM[z-base,:]=[z,np.mean(location[1][findz[0]]), np.mean(location[2][findz[0]])] #find center of mass w/ arithmetic mean
            
        
        area_z_ratio[m]= np.mean(cross_area_z)/np.amax(cross_area_z)    
        center_choice=[np.mean(COM[:,1]), np.mean(COM[:,2])] ### choosing pt to align COMs
        movement=np.subtract(center_choice , COM[:,1:])
        for z in np.arange(base,top+1,1): ### shifting cld
            findz=np.where(location[0] == z)
            newlocation_matrix[findz,1:]=location_matrix[findz,1:] + movement[z-base,:]
        
        newlocation_matrix[:,0]=location_matrix[:,0]
        
        location[0][:]=newlocation_matrix[:,0]   # change to calculate overlap w/ code below
        location[1][:]=newlocation_matrix[:,1]
        location[2][:]=newlocation_matrix[:,2]
        
        #####################################\
        ### calculate shift distance b/c => w = (wind diff * ht)/ shift_distance
        movement1 = 25*COM[:,1:] ### having center at origin,spacing in x and y direction
        shift_cld= np.sqrt(movement1[:,0]**2 + movement1[:,1]**2)
        #shift_cld = (abs(movement[:,0])  + abs(movement[:,1]) ) / 2
        shift_max[m] = np.amax(shift_cld)
        shift_min[m] = np.amin(shift_cld)
        shift_avg[m] = np.mean(shift_cld)
        shift_distance[m] = np.amax(shift_cld) - np.amin(shift_cld) 
        #########################################
        
        ########
        
        height_c[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
        #projected_area_c[m]=dx*dy*len(np.unique(location[2]+location[1]*nx)) #len of unique # of cells that get proj
        partial_heights=alpha*height_c[m] # fraction of cld height
        #height_diff=height_c[m]-partial_heights
        
        overlap_changed=np.zeros(partial_heights.size)
        projected_area_changed=np.zeros(partial_heights.size)
        vol_changed=np.zeros(partial_heights.size)
        part_loc= -1
        for partial in partial_heights: ### probably dont need this loop since we are not chopping cld ht
            
            part_loc= part_loc +1 # to index later for diff ht levels
            for i in range(location[0].size):
                if dz*(location[0][i]-location[0][0]+1) <= partial: # finding all cells below chopped cld ht (take highest cell)
                    celltop=i
                else:
                    break
            #if partial - dz*(location[0][celltop]-location[0][0]+1) >0: 
            difference = partial - dz*(location[0][celltop]-location[0][0]+1) # diff btw chopped cld ht and ht of clds w/ full cells
            #print(difference == 0 )
            #print(difference)
            for j in range(celltop,location[0].size):
                if  location[0][celltop] +1  == location[0][j]: # to use partial cells
                    cellabove=j
                elif location[0][celltop] +1  < location[0][j]:
                    break
            #print('partial',partial)    
            #print('top',celltop)
            #print('above',cellabove)
            #loc_loc0=-1
            if partial == height_c[m]: # difference=0, so stop at celltop
                proj_area_loop1=np.zeros(celltop +1)
                vol_loop1=np.zeros(celltop +1)
                
                """
                ##### 2d plot of cld
                bins=dz
                plt.figure()
                plt.hist2d(dx*location[2][0:celltop+1],dz*location[0][0:celltop+1],bins=bins,cmin=0.5)
                plt.xlabel('x');plt.ylabel('z')
                plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                colorbar = plt.colorbar()
                colorbar.set_label('counts in bin')
                plt.figure()
                plt.hist2d(dy*location[1][0:celltop+1],dz*location[0][0:celltop+1],bins=bins,cmin=0.5)
                plt.xlabel('y');plt.ylabel('z')
                plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                colorbar = plt.colorbar()
                colorbar.set_label('counts in bin')
                """
                """
                 ### plotting in cld in 3d
        
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.scatter3D(dx*location[2][0:celltop+1],dy*location[1][0:celltop+1], dz*location[0][0:celltop+1], c=dz*location[0][0:celltop+1] ,cmap='Paired'); 
                #cmap= single color:copper, cool, winter, multicolor:  Dark2, Paired
                #ax.plot_trisurf(location[2], location[1], location[0], cmap='viridis', edgecolor='none');
                plt.title(str(cloud_numb[m]))
                plt.show()
                """
                
            else:
                proj_area_loop1=np.zeros(cellabove +1)
                vol_loop1=np.zeros(cellabove +1)
                
                """
                ##### 2d plot of cld
                bins=dz
                plt.figure()
                plt.hist2d(dx*location[2][0:cellabove+1],dz*location[0][0:cellabove+1],bins=bins,cmin=0.5)
                plt.xlabel('x');plt.ylabel('z')
                plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                colorbar = plt.colorbar()
                colorbar.set_label('counts in bin')
                plt.figure()
                plt.hist2d(dy*location[1][0:cellabove+1],dz*location[0][0:cellabove+1],bins=bins,cmin=0.5)
                plt.xlabel('y');plt.ylabel('z')
                plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                colorbar = plt.colorbar()
                colorbar.set_label('counts in bin')
                """
                        
                
            for loc0 in range(location[0].size): ### calculating proj area and vol
                
                #loc_loc0 = loc_loc0 +1
                if loc0  <= celltop : 
                    proj_area_loop1[loc0]=(location[2][loc0]+location[1][loc0]*nx)
                elif celltop < loc0 and loc0 <=cellabove:
                    proj_area_loop1[loc0]=(location[2][loc0]+location[1][loc0]*nx)
                if loc0  <= celltop: # full cell vol
                    vol_loop1[loc0]=(dx*dy*dz)   
                elif celltop < loc0 and loc0 <=cellabove: # partial cell vol
                    vol_loop1[loc0]=(dx*dy*difference) 
                else:
                    break
                #if dz*loc0 < parital:
            #proj_area_loop1=np.asarray(proj_area_loop1)
            #print('loop1',proj_area_loop1.size)
            projected_area_changed[part_loc]=dx*dy*len(np.unique(proj_area_loop1)) # counting cells w/ ! id #  & proj-ing
            #vol_loop1=np.asarray(vol_loop1)  
            vol_changed[part_loc]=sum(vol_loop1)
            overlap_changed[part_loc]=vol_changed[part_loc]/(partial_heights[part_loc]*projected_area_changed[part_loc])
            if partial == height_c[m]:
                overlap_c[m]=overlap_changed[part_loc] # make sure calculated overlap is correct for full ht
                volume_c[m]=vol_changed[part_loc]
                projected_area_c[m]=projected_area_changed[part_loc]
                
            
        
        
        #m = m+1
    
    #mu1=np.mean(Total_dr)
    #sigma1=np.std(Total_dr)
    #Check1=[Total_dr >= mu1 + 1.96 *sigma1]
    
    overlap_c=volume_c/(height_c*projected_area_c)
    
    
    
    
    """
    ####################################################
    ##### saving data for convex hull
    #if dz*layers >100:
    #    overlap_contra=1/overlap_ratio - 1/overlap_convex 
        
    
    #step=1800 # 30*60
    #time_array= np.arange(step,36000+step,step)
    #npyfilespath ="/home/anthonys/Documents/bomex_overlap_convex"
    ###
    npyfilespath ="/home/anthonys/Documents/"
    
    name1 = 'overlap_convex_'
    
    #np.save(npyfilespath+str(time_array[file1_numb])+'.npy',overlap_convex)
    np.save(npyfilespath+ name1  + case +str(time_array[file1_numb])+'.npy',overlap_convex)
    
    ############################################################
    
    #### save wind differentials and updraft vel.
    
    
    #step=1800 # 30*60
    #time_array= np.arange(10800,36000+step,step)
    npyfilespath ="/home/anthonys/Documents/"
    #case='bomex'
    
    
    name1='shearTB_'
    name2='shear_sum_'
    #name3='wavg1_'
    
    
    np.save(npyfilespath+ name1  + case +str(time_array[file1_numb])+'.npy',shear_cal)
    np.save(npyfilespath+ name2  + case +str(time_array[file1_numb])+'.npy',shear00)
    #np.save(npyfilespath+ name3  + case +str(time_array[file1_numb])+'.npy',w_avg)
    
    name1='wz_max_'
    name2='wz_min_'
    name3='wz_avg_'
    name4='wz_ct_'
    name5='wz_cb_'
    
    
    np.save(npyfilespath+ name1  + case +str(time_array[file1_numb])+'.npy',wz_max)
    np.save(npyfilespath+ name2  + case +str(time_array[file1_numb])+'.npy',wz_min)
    np.save(npyfilespath+ name3  + case +str(time_array[file1_numb])+'.npy',wz_avg)
    np.save(npyfilespath+ name4  + case +str(time_array[file1_numb])+'.npy',wz_ct)
    np.save(npyfilespath+ name5  + case +str(time_array[file1_numb])+'.npy',wz_cb)
    
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
    
    
    
    ################################
    ### saving data 
    


    #step=1800 # 30*60
    #time_array= np.arange(step,36000+step,step)
    npyfilespath ="/home/anthonys/Documents/"
    #case='bomex'
    ###
    #step=3600 # 60*60
    #time_array= np.arange(step,216000+step,step)
    #npyfilespath ="/home/anthonys/Documents/"
    #case='rico'
    ###
    #step = 600
    #time_array = np.zeros(56)
    #time_array[0:2] = [0,1200,6600]
    #time_array[3:] = np.arange(7800,41400+step,step)
    npyfilespath ="/home/anthonys/Documents/"
    #case='lasso'
    
    name1='overlap_no_shear_'
    name2='volume_no_shear_'
    name3='projarea_no_shear_'
    name4='height_no_shear_'
    name5='area_z_ratio_'
    
    
    ### 
    np.save(npyfilespath+ name1  + case +str(time_array[file1_numb])+'.npy',overlap_c)
    np.save(npyfilespath+ name2  + case +str(time_array[file1_numb])+'.npy',volume_c)
    np.save(npyfilespath+ name3  + case +str(time_array[file1_numb])+'.npy',projected_area_c)
    np.save(npyfilespath+ name4  + case +str(time_array[file1_numb])+'.npy',height_c)
    np.save(npyfilespath+ name5  + case +str(time_array[file1_numb])+'.npy',area_z_ratio)
    
    
    #np.save(npyfilespath+str(time_array[20])+'.npy',area_z_ratio)

    

    ###########################
    """













###################################

end= time.time()
print('Run Time in Seconds:', end-start)

