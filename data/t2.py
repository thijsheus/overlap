#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 10:59:43 2019

@author: anthonys
"""






#from netCDF4 import Dataset
import numpy as np 
#import struct 
#import netCDF4
from netCDF4 import Dataset
#import collections
import matplotlib.pyplot as plt
#from scipy.io import netcdf
#import scipy as sp
#import itertools
import os
import sys
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit



#############################################################################
bomexd = Dataset("bomex.default.0000000.nc","r")
armd= Dataset("arm.default.0000000.nc","r")
ricod= Dataset("rico.default.0000000.nc","r")


#print(bomexd.variables.keys())
#print(bomexd.variables['area'])
#print(bomexd.variables['area'][:][:])
#print(bomexd.variables)
#print(bomexd.variables.values())
#########################################################################
# Define variables and importing files

area  = bomexd.variables["area"][:,:] # fractional area contained in mask
areah = bomexd.variables['areah'][:,:]
time=bomexd.variables['time'][:] 

z=bomexd.variables['z'][:] # full level ht 
zh=bomexd.variables['zh'][:] # half level height
ql=bomexd.variables['ql'][:,:] 

thl=bomexd.variables['thl'][:,:] # liquid water potential temp
qlfrac=bomexd.variables['qlfrac'][:,:] # cloud fraction
qlcover=bomexd.variables['qlcover'][:] # projected cloud cover

u2=bomexd.variables['u2'][:,:]  #moment 2 of u vel, weighted ave of intenties of pixels in seq of images
u3=bomexd.variables['u3'][:,:]
u4=bomexd.variables['u4'][:,:]

ugrad=bomexd.variables['ugrad'][:,:] # gradient (shear) of the u velocity
u=bomexd.variables['u'][:,:] # vel in u direction
v=bomexd.variables['v'][:,:] # vel in v direction
w=bomexd.variables['w'][:,:] # vertical velocities



b=bomexd.variables['b'][:,:] # buoyancy
rho=bomexd.variables['rho'][:,:] # full level density
rhoh=bomexd.variables['rhoh'][:,:] # half level density
phydro=bomexd.variables['phydro'][:,:] # full level hydrostatic pressure

"""
bomexql=Dataset("bomex.ql.0000000.nc","r")
qlfrac1=bomexql.variables['qlfrac'][:,:]
z1=bomexql.variables['z'][:,:]
zh1=bomexql.variables['zh'][:,:]


bomexqlcore=Dataset("bomex.qlcore.0000000.nc","r")
qlfrac2=bomexqlcore.variables['qlfrac'][:,:]
z2=bomexqlcore.variables['z'][:]
zh2=bomexqlcore.variables['zh'][:]
"""

bomextrack=Dataset("l.0001800.track.nc","r") 
cv=bomextrack.variables['cv'][:] # cloud cover by volume 
cp=bomextrack.variables['cp'][:] # cloud cover by area
area_proj=bomextrack.variables['area_proj'][:] # projected cloud area
cloudsize=np.sqrt(area_proj)
area_proj_conn=bomextrack.variables['area_proj_conn'][:] # connected projected cloud area
ht=bomextrack.variables['ht'][:] # height
z=bomextrack.variables['z'][:]
x=bomextrack.variables['x'][:]
y=bomextrack.variables['y'][:]
cld_mask=bomextrack.variables['cld_mask'][:,:,:] # 
cfrac=bomextrack.variables['cfrac'][:] # cloud fraction:percentage of each pixel in gridbox in model that is covered w/ clouds
nrcloud=bomextrack.variables['nrcloud'][:,:,:] # cloud number, masks (inte array)
nr=bomextrack.variables['nr'][:] # number density rain ?
cb=bomextrack.variables['cb'][:] # cloud base
ct=bomextrack.variables['ct'][:] # cloud top,  # ct - cb , height difference
overlap_ratio=bomextrack.variables['chr'][:] # overlap ratio based on cloud height 
cfv=bomextrack.variables['cfv'][:,:] # cv for cloud layer/field, cloud field projected area 
cfp=bomextrack.variables['cfp'][:,:] # cp for cloud layer/field, cloud field vol per ht
# cfr=cfv/cfp # overall ratio
# rav=bomextrack.variables['rav'][:] # average ratio

nz=z.size;nx=x.size;ny=y.size; # sizes of dimensions




bomextrack1=Dataset("l.0003600.track.nc","r") 
"""
cv1=bomextrack1.variables['cv'][:]
cp1=bomextrack1.variables['cp'][:]
"""

bomextrack2=Dataset("l.0005400.track.nc","r") 
"""
cv2=bomextrack.variables['cv'][:]
cp2=bomextrack.variables['cp'][:]
area_proj2=bomextrack.variables['area_proj'][:]
cloudsize2=np.sqrt(area_proj)
ht2=bomextrack.variables['ht'][:]
"""

bomdaleprof=Dataset("profiles.001.nc","r")


bomexdalesamp=Dataset('sampling.001.nc','r')



###########################################################################


ht_sort=np.sort(ht)
ht_sort=np.float64(ht_sort) ### needs to match type (e.g. float64) w/ codomain
#ht.argsort()
overlap_ratio_sort=np.zeros(int(overlap_ratio.size))
j=0
for i in ht.argsort():
    overlap_ratio_sort[j]=overlap_ratio[i]
    j=j+1



def invlin(x,a,b):
    return 1/(a+b*x)


def plaw(x,a,b):
    return a*x**b


vertaxis=ht_sort  # defining dummy var to easily change var 
horaxis=1/overlap_ratio_sort    
parm, parmcov=curve_fit(plaw,vertaxis,horaxis)
fit=plaw(vertaxis,*parm)  

"""
plt.figure()
plt.plot(vertaxis,overlap_ratio_sort,'o')
plt.plot(vertaxis,fit,'-')
"""
plt.figure()
#plt.plot(overlap_ratio_sort,vertaxis,'o')
plt.hist2d(horaxis, vertaxis, bins=15,cmin=0.5)
cb = plt.colorbar()
cb.set_label('counts in bin')
plt.plot(fit,vertaxis,'-')
plt.title('Power Law: ht vs. 1/overlap',fontsize=18)
print('Power Law parameters:',parm)


#################################################

ht_sort=np.sort(area_proj)
ht_sort=np.float64(ht_sort) ### needs to match type (e.g. float64) w/ codomain




vertaxis=ht_sort  # defining dummy var to easily change var 
horaxis=1/overlap_ratio_sort   
parm, parmcov=curve_fit(plaw,vertaxis,horaxis)
fit=plaw(vertaxis,*parm)  


plt.figure()
#plt.plot(overlap_ratio_sort,vertaxis,'o')
plt.hist2d(horaxis, vertaxis, bins=15,cmin=0.5)
cb = plt.colorbar()
cb.set_label('counts in bin')
plt.plot(fit,vertaxis,'-')
plt.title('Power Law: area_proj vs. 1/overlap',fontsize=18)
print('Power Law parameters:',parm)



































