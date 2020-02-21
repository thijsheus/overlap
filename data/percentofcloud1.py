#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 11:56:08 2019

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


filenames=[bomexd, ricod, armd]

bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72, bomextrack90, bomextrack108, bomextrack126, bomextrack144, bomextrack162, bomextrack180, bomextrack198, bomextrack216]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144, ricotrack180, ricotrack216, ricotrack252, ricotrack288, ricotrack324, ricotrack360, ricotrack396]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162, armtrack180, armtrack198, armtrack216, armtrack234, armtrack252, armtrack270, armtrack288]

###########################################################################



# definitions
def pyth(u,v):    # magnitude
    return np.sqrt(u*u+v*v)

########
# when |ct-cb|<1 ra -> rashear
def rashear(grad,ht,w,ct,cb,one):
    multiple=(abs(grad)*ht)/np.maximum(abs(w),one)
    multiplicand=1/cb
    return (1+multiple*multiplicand)**(-1)

def raanvil(grad,ht,w,ct,cb,one):
    multiple=1/np.maximum(abs(w),one)
    multiplicand=np.maximum(abs(ct-cb),one)/cb
    return (1+multiple*multiplicand)**(-1)

# model
def ra(grad,ht,w,ct,cb,one):
    multiple=(np.maximum(abs(grad)*ht,one))/np.maximum(abs(w),one)
    multiplicand=np.maximum(abs(ct-cb),one)/cb
    return ((1)+multiple*multiplicand)**(-1)





####################################
    
# script

bomexfilenames=[bomextrack18]



for file in filenames:
    #zt=file.variables['z'][:]
    zh=file.variables['zh'][:]
    time_t=file.variables['time'][:]
    ugrad=file.variables['ugrad'][:,:]
    vgrad=file.variables['vgrad'][:,:]
    u=file.variables['u'][:,:]
    v=file.variables['v'][:,:]
    w=file.variables['w'][:,:]

    
            
    if file == bomexd:
        for file1 in bomexfilenames:
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
            
            
            
            cbnew1=[];ctnew1=[];cbnew2=[];ctnew2=[];area_projnew1=[];area_projnew2=[];
            overlap_ratio_new1=[];overlap_ratio_new2=[];htnew1=[];htnew2=[];
            for i in range(cb.size):
                if abs(ct[i]-cb[i]+1) > 5:
                    cbnew1.append(cb[i])
                    ctnew1.append(ct[i])
                    area_projnew1.append(area_proj[i])
                    overlap_ratio_new1.append(overlap_ratio[i])
                    htnew1.append(ht[i])
                else:
                    cbnew2.append(cb[i])
                    ctnew2.append(ct[i])
                    area_projnew2.append(area_proj[i])
                    overlap_ratio_new2.append(overlap_ratio[i])
                    htnew2.append(ht[i])
                
            
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

            
            
            
            
            
            
            
            
            
            dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
            gridarea=dx*dy
            gridvol=dx*dy*dz
            bins=5;
            """
            for z in range(zt.size):
                hist1,edges = np.histogram(nrcloud[z,:,:],bins=bins)
            """

            
            nrcloudarray = np.ma.getdata(nrcloud)
            cld_maskarray = np.ma.getdata(cld_mask)
            """
            cloudfrac = np.empty((zt.size))
            for z in range(zt.size):
                cloudfrac[z] = np.sum((nrclouddata[z,:,:]))/(xt.size*yt.size)
            """
            
            """
            cvlayer=np.zeros(cfrac.size)
            for i in range(cfrac.size):
                cvlayer[i]=cfrac[i]/(dz+1)
            """
            
            """
            N=[]
            for i in range(xt.size):
                for j in range(yt.size):
                    #Z=sum(nrcloudarray[:,j,i])
                    #N.append(Z)
                    for k in range(zt.size):
                        if nrcloudarray[k,j,i]>0:
                            N.append(1)
                    
            Narray=np.asarray(N)
            S=sum(Narray)
            """

            
            count1=np.count_nonzero(nrcloudarray)
            count2=np.count_nonzero(cld_maskarray)

            #nrcloudarray.index(1)
            #nrcloudlist=np.tolist(nrcloudarray)
            #nrcloudarray.index(1)
            A=np.where(nrcloudarray==158)


            
















end= time.time()
print('Run Time in Seconds:', end-start)
















