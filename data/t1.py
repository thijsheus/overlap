#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:21:13 2019

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



# plotting


"""
plt.figure(1)
plt.plot(ql,thl)
plt.xlabel('Liquid water mixing ratio')
plt.ylabel('Liquid water potential temperature')




plt.figure(2)
plt.plot(qlfrac)
plt.title('Cloud Fraction')


plt.figure(3)
plt.plot(thl,qlfrac,'o')
plt.title('qlfrac vs. thl')


plt.figure(4)
[rows,columns]=np.shape(qlfrac)
for i in range(rows):
    plt.plot(z,qlfrac[i,:])
plt.title('cloud fraction vs. z ')



plt.figure(5)
plt.plot(u2,u3,u4,'o')
plt.title('moments in the u direction')


plt.figure(6)
plt.plot(u,v,w[:,0:-1],'o')
plt.title('Velocities')
plt.ylim([-0.05, 0.15])



plt.figure(7)
[rows,columns]=np.shape(b)
for i in range(rows):
    plt.plot(z,b[i,:])
plt.title('buoyancy vs. z ')

plt.figure(8)
plt.plot(time,qlcover);
plt.title('cloud cover vs. time')



plt.figure(9)
[rows,columns]=np.shape(rho)
for i in range(rows):
    plt.plot(z,rho[i,:])
plt.title('Full Level Densities vs. z ')
  


plt.figure(10)
plt.plot(rho,phydro)
plt.title('Full Level Densities vs. Full Level Hydrostatic Pressure')


plt.figure(11)
plt.plot(zh[0:-1],z,'*')
plt.title('Full level height vs. half level height')
"""





# Plots with Overlap Ratio

"""
plt.figure()
plt.plot(overlap_ratio,area_proj,'o')

plt.title('Cloud Overlap Ratio vs. Area', fontsize=18)


plt.figure()
plt.plot(overlap_ratio,cloudsize,'o')
plt.title('Cloud Overlap Ratio vs. Cloud Size', fontsize=18)


plt.figure()
plt.plot(overlap_ratio,ht,'o')
plt.title('Cloud Overlap Ratio vs. ht', fontsize=18)
"""

# 3D plot
"""
ax = plt.axes(projection='3d')
ax.plot3D(ht,area_proj,overlap_ratio,'o')
ax.set_xlabel('ht', fontsize=14)
ax.set_ylabel('area', fontsize=14)
ax.set_zlabel('cv/cp', fontsize=14)
ax.set_title('cv/cp vs. ht & area',fontsize=18)
"""

# Curve fit proj area vs. overlap ratio w/ inverselinear 
"""
def invlin(x,b):
    return 1/(1+b*x)

vertaxis=ht  # defining dummy var to easily change var 

parm,parmcov=curve_fit(invlin,vertaxis,overlap_ratio)
fit=invlin(vertaxis,*parm)  # need monotonic domain
#plt.figure()
#plt.plot(area_proj,overlap_ratio,'o')
#plt.plot(area_proj,invlin(area_proj,*parm),'*')
plt.figure()
plt.plot(overlap_ratio,vertaxis,'o')
plt.plot(fit,vertaxis,'-')
"""


"""
ugrad=bomexd.variables['ugrad'][:,:] # gradient (shear) of the u velocity
w=bomexd.variables['w'][:,:] # vertical velocities

##h=ht;b=cb;u=ugrad;
#r_b=1/(1+(ugrad/(2*w))+(ht/cb)) # need to convert to same size
"""
"""
zdiff=np.zeros(int(z.size-1)) # z is evenly spaced
for i in range(int(z.size-1)):
    zdiff[i]=z[i+1]-z[i]
"""

 #plt.plot(z,cfrac,'*') 


"""
plt.figure()
plt.plot(cloudsize,ht,'o')
plt.title('cloud size vs. ht', fontsize=18)
plt.figure()
plt.plot(z,cfrac,'o')
plt.title('z vs. cfrac', fontsize=18)
"""

################################################################################
# sorting input to get strictly monotonic domain and curve fitting 

ht_sort=np.sort(ht)
ht_sort=np.float64(ht_sort) ### needs to match type (e.g. float64) w/ codomain
#ht.argsort()
overlap_ratio_sort=np.zeros(int(overlap_ratio.size))
j=0
for i in ht.argsort():
    overlap_ratio_sort[j]=overlap_ratio[i]
    j=j+1

#plt.figure()
#plt.plot(overlap_ratio_sort,ht_sort,'o')


def invlin(x,a,b):
    return 1/(a+b*x)


def plaw(x,a,b):
    return a*x**b


vertaxis=ht_sort  # defining dummy var to easily change var 
    
parm, parmcov=curve_fit(plaw,vertaxis,overlap_ratio_sort)
fit=plaw(vertaxis,*parm)  

"""
plt.figure()
plt.plot(vertaxis,overlap_ratio_sort,'o')
plt.plot(vertaxis,fit,'-')
"""
plt.figure()
plt.plot(overlap_ratio_sort,vertaxis,'o')
plt.plot(fit,vertaxis,'-')
plt.title('Power Law: ht vs. overlap',fontsize=18)
print('Power Law parameters:',parm)


parm, parmcov=curve_fit(invlin,vertaxis,overlap_ratio_sort)
fit=invlin(vertaxis,*parm)  

plt.figure()
plt.plot(overlap_ratio_sort,vertaxis,'o')
plt.plot(fit,vertaxis,'-')
plt.title('Inverse Linear: ht vs. overlap',fontsize=18)
print('Inverse Linear parameters:',parm)

##############################################################################

#######################################################################
 #sorting input to get strictly monotonic domain and curve fitting 

ht_sort=np.sort(area_proj)
ht_sort=np.float64(ht_sort) ### needs to match type (e.g. float64) w/ codomain
#ht.argsort()
overlap_ratio_sort=np.zeros(int(overlap_ratio.size))
j=0
for i in area_proj.argsort():
    overlap_ratio_sort[j]=overlap_ratio[i]
    j=j+1

#plt.figure()
#plt.plot(overlap_ratio_sort,ht_sort,'o')


def invlin(x,a,b):
    return 1/(a+b*x)


def plaw(x,a,b):
    return a*x**b


vertaxis=ht_sort  # defining dummy var to easily change var 
    
parm, parmcov=curve_fit(plaw,vertaxis,overlap_ratio_sort)
fit=plaw(vertaxis,*parm)  


plt.figure()
plt.plot(overlap_ratio_sort,vertaxis,'o')
plt.plot(fit,vertaxis,'-')
plt.title('Power Law: area_proj vs. overlap',fontsize=18)
print('Power Law parameters:',parm)

"""
parm, parmcov=curve_fit(invlin,vertaxis,overlap_ratio_sort)
fit=invlin(vertaxis,*parm)  

plt.figure()
plt.plot(overlap_ratio_sort,vertaxis,'o')
plt.plot(fit,vertaxis,'-')
plt.title('Inverse Linear: area_proj vs. overlap',fontsize=18)
print('Inverse Linear parameters:',parm)
"""












#####################################################################








          
##############################################################################
plt.show()



















