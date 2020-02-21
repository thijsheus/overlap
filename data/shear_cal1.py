#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:55:59 2019

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
import matplotlib.colors as colors
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import hmean
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

from scipy.special import comb


plt.rcParams.update({'font.size': 16})

start=time.time()


#############################################################################

# definitions
#def pyth(u,v):    # magnitude
#    return np.sqrt(u*u+v*v)
    


#### definitions
def f_turb(h,parm):
    return (3*np.pi /4)* (h/(parm+h))

def f_shear(aspect,parm):
    return parm*aspect**1

def f_area(aspect,parm):
    return parm*aspect**1

def shear_amod(shear,h,l):
    return shear * h / l


####################################
    
# script

epsilon1=50    #arbitrarily set, most outliers are > 200, most nonoutliers clds are about 0

#kt=0;
#file_numb=-1
#file1_numb=-1



### rico828
"""
height=np.load('height_no_shear_rico828.npy')
volume=np.load('volume_no_shear_rico828.npy')
projar=np.load('projarea_no_shear_rico828.npy')
overlap=np.load('overlap_no_shear_rico828.npy')
areaz=np.load('area_z_ratio_rico828.npy')
areazg=np.load('rico_area_gmean82800.npy')
overlap_convex=np.load('overlap_convex_rico82800.npy')
"""
"""
### rico2016
height=np.load('height_no_shear_rico2016.npy')
volume=np.load('volume_no_shear_rico2016.npy')
projar=np.load('projarea_no_shear_rico2016.npy')
overlap=np.load('overlap_no_shear_rico2016.npy')
areaz=np.load('area_z_ratio_rico2016.npy')
model=np.load('model_rico2016.npy')
overlap_convex=np.load('overlap_convex_rico201600.npy')
"""

"""
### bomex 360
height=np.load('height_no_shear_bomex360.npy')
volume=np.load('volume_no_shear_bomex360.npy')
projar=np.load('projarea_no_shear_bomex360.npy')
overlap=np.load('overlap_no_shear_bomex360.npy')
areaz=np.load('area_z_ratio_bomex360.npy')
"""
"""
### lasso 306
height=np.load('height_no_shear_lasso306.npy')
volume=np.load('volume_no_shear_lasso306.npy')
projar=np.load('projarea_no_shear_lasso306.npy')
overlap=np.load('overlap_no_shear_lasso306.npy')
areaz=np.load('area_z_ratio_lasso306.npy')
"""


"""
datad = Dataset("/data/rico/rico.default.0000000.nc","r")
zh=datad.variables['zh'][:]
time_t=datad.variables['time'][:]
u=datad.variables['u'][:,:]
v=datad.variables['v'][:,:]
w=datad.variables['w'][:,:]
"""
####################################################################
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

for file1 in Bfilenames:
    print(file1)

    data = Dataset(file1,'r')

    #file1_numb = file1_numb+1
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
    
    
    ###
    access = np.where(time_t == time_array[file1_numb])[0][0] ## slicing at correct spot
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
    
    """
    ### Bomex 
    overlap_ratio=np.load('Bomex_overlap_ratio.npy')            
    height=np.load('Bomex_ht.npy')
    overlap=np.load('Bomex_shear.npy')
    areaz=np.load('Bomex_areaz.npy')
    ht=height
    cv=np.load('Bomex_cv.npy')
    cp=np.load('Bomex_cp.npy')
    areazg=np.load('Bomex_area_gmean.npy')
    overlap_convex=np.load('Bomex_hull.npy')
    print('Bomex')
    """
    
    l=cv**0.5
    
    
    
    
    
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
    
    
    
    ##############################
    
    duz=np.gradient(uz) # central diff in uz
    dvz=np.gradient(vz)
    
    
    
    #s_u=abs(duz);s_v=abs(dvz)
    s_u=duz;s_v=dvz;
    
    #Gmean=(s_u*s_v)**(0.5)
    Amean=(s_u+s_v)/2;#Amean=abs(Amean)
    #H=Gmean**2/Amean
    Hmean=(2*s_u*s_v) / (s_u + s_v);#Hmean=abs(Hmean)
    pyth= (s_u**2 + s_v**2)**0.5
    
    shear0=np.zeros(cb.size)
    shear00=np.zeros(cb.size)
    w_cld=np.zeros(cb.size)
    for i in range(cb.size): #### wind differential 
        
        #shear0[i] = sum(Hmean[cb[i]:ct[i]+1])    # since max(w,1)=1 s_u,s_v = duz, dvz
        shear00[i] = sum(Amean[cb[i]:ct[i]+1])
        
        shear_calu[i] = (uz[ct[i]] - uz[cb[i]])
        shear_calv[i]= (vz[ct[i]] - vz[cb[i]])
        
    shear_cal = abs(shear_calu + shear_calv) / 2
    
    #######################
    #shear=   (1/overlap_ratio) - (1/overlap) ### overlap form shear
    
    index1=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
    index1_size=index1[0].size
    ht=ht[index1[0]]   # taking the ht values according to indices above
    overlap_ratio=overlap_ratio[index1[0]];
    area_proj=area_proj[index1[0]]
    print('Clouds that satisfy the given condition: ',index1_size)
    print('conditional height is',conditional_height,'meters')

    ll = np.zeros(index1_size)
    shear1_cal = np.zeros(index1_size)
    shear2_cal = np.zeros(index1_size)
    shear3_cal = np.zeros(index1_size)
    w_avg = np.zeros(index1_size)
    shear_overlap=np.zeros(index1_size)
    
    C_true_bar = np.zeros(index1_size)
    alpha_overlap = np.zeros(index1_size)
    cloud_numb=index1[0] +1 # index is off by 1 as array starts with zero
    
    
    wz_max = np.zeros(index1_size)
    wz_cb = np.zeros(index1_size)
    wz_ct = np.zeros(index1_size)
    wz_avg = np.zeros(index1_size)
    wz_min = np.zeros(index1_size)
    m=-1; #m=0;
    
    
    #### calculating updraft
    for e1 in cloud_numb:  # may take a while to run 
        
        m = m+1
        
        ll[m]=l[e1-1]
        
        location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1
    
        location_matrix=np.zeros((location[0].size,3))
        location_matrix[:,0]=location[0] # z coor
        location_matrix[:,1]=location[1] # y coor
        location_matrix[:,2]=location[2] # x coor
        
        loc_mat_tran=np.matrix.transpose(location_matrix)
        layers=(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z 
        base=np.amin(location[0]);top=np.amax(location[0]);
        z_array = np.arange(base, top+1,1)
        ###########
        ### shear
        basez=base;topz=top;basey=min(location[1]);topy=max(location[1]);basex=min(location[2]);topx=max(location[2]);
        w_abs=abs(wt)
        #wz_avg=np.sum(w_abs[basez:topz+1,:,:],axis=0) / (topz - basez+1)  # maybe avoid avg along z axis? use max up vel along z to get singl #
        #wy_avg=np.sum(wz_avg[basey:topy+1,:],axis=0) / (topy - basey+1)
        #wx_avg=np.sum(wy_avg[basex:topx+1],axis=0) / (topx - basex+1)
        ##wz_avg=np.mean(w_abs[min(location[0]):max(location[0])+1,:,:],axis=0)
        ##wy_avg=np.mean(wz_avg[min(location[1]):max(location[1])+1,:],axis=0)
        ##wx_avg=np.mean(wy_avg[min(location[2]):max(location[2])+1],axis=0)
        
        wy_avg=np.sum(w_abs[:,basey:topy+1,:],axis=1) / (topy - basey+1)
        wx_avg=np.sum(wy_avg[:,basex:topx+1],axis=1) / (topx - basex+1)
        wz = wx_avg[basez:topz+1]
        
        wz_max[m] = np.amax(wz)
        wz_cb[m] = wz[0]
        wz_ct[m] = wz[-1]
        wz_avg[m] = np.mean(wz)
        wz_min[m] = np.amin(wz)
        
        
        #w_avg[m]=wx_avg
        #print(w_avg[m])
        
        
        
        #shear1_cal[m]= shear_cal[e1-1]/ w_avg[m] ###shear_cal
        
        #shear2_cal[m] = shear0[e1-1] / w_avg[m] ### shearh
        
        #shear3_cal[m] = shear00[e1-1] / w_avg[m] ###sheara
        #print(shear1_cal)
    
    
        #shear_overlap[m] = shear[e1-1]
        #########################
        ### hogan ,illington
        """
        area_layer = np.zeros(z_array.size)
        for z in z_array:
            area_layer[z-base] = dx*dy*len(np.where( z == location[0] )[0])
        
        
        C_max = np.zeros(int(comb(z_array.size,2)))
        C_rand = np.zeros(int(comb(z_array.size,2)))
        C_true = np.zeros(int(comb(z_array.size,2)))
        alpha_a = np.zeros(int(comb(z_array.size,2)))
        loc_a_1 = -1
        for c_area1 in area_layer:
            loc_a_1 = loc_a_1 + 1
            loc_a_2 = - 1
            for c_area2 in area_layer:
                loc_a_2 = loc_a_2 + 1
                if loc_a_1 != loc_a_2:
                    dzed = dz*abs(loc_a_1 - loc_a_2)
                    C_max[loc_a_2] = max(c_area1,c_area2)
                    C_rand[loc_a_2] = c_area1 + c_area2 - c_area1*c_area2 
                    alpha_a[loc_a_2] = np.exp(-dzed / 200 )
                    C_true[loc_a_2] = alpha_a[loc_a_2]*C_max[loc_a_2] + (1-alpha_a[loc_a_2])*C_rand[loc_a_2]
                    
        #C_max_bar = np.mean(C_max)
        #C_rand_bar = np.mean(C_rand)
        C_true_bar[m] = np.mean(C_true)
        alpha_overlap[m] = np.exp(-dz*layers / 200)
        """
        
    ###################################################
    """
    fract= (5*np.pi) / (6*np.sqrt(3))  -1
    fractal=  (5*np.pi) / (6*np.sqrt(3)) # -1
    fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) -1
    
    shape=  (1/areaz) - 1
    shear=   (1/overlap_ratio) - (1/overlap)  
    shear[shear<0]=0
    
    turb= (1/overlap_ratio[ht>100]) - (1/overlap_convex)
    turb[turb<0]=0
    fractalf[ht>100]=turb
    
    
    ## upper bound to turb
    lam =   (3*np.pi)/4  #np.pi * np.sqrt(3) / 2     #(15/8)*np.sqrt(3/8)*np.pi
    #C1=np.argwhere((1/overlap_ratio)>1.5*lam)
    keep=np.argwhere(fractalf > lam)
    fractalf[np.array(keep)] = lam
    
    
    r_total =  ( (shear) + (shape)  + (fractalf) ) +1
    
    
    f_t=f_turb(ht,175); f_t[ht<=100]=0
    f_s=f_shear(ht/l,0.43);f_s[ht/l<=1]=0
    f_a=f_area(ht/l,0.2);f_a[ht/l<=1]=0
    
    #r_total = f_t + f_s + f_a +1
    
    ##################################################
    """
    ################################
    ### saving data
    """
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
    """
    #############################
    
###############################
                
end= time.time()
print('Run Time in Seconds:', end-start)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    