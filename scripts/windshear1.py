#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 16:06:23 2019

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


start=time.time()

#############################################################################

armd= Dataset("arm.default.0000000.nc","r")
armql=Dataset("arm.ql.0000000.nc","r")
armqlcore=Dataset("arm.qlcore.0000000.nc","r")
ricod= Dataset("rico.default.0000000.nc","r")
ricoql=Dataset("rico.ql.0000000.nc","r")
ricoqlcore=Dataset("rico.qlcore.0000000.nc","r")

#data = pkgutil.get_data("conf", "fielddump.ql.00.05.track.nc")


#print(bomexd.variables.keys())
#print(bomexd.variables['area'])
#print(bomexd.variables['area'][:][:])
#print(bomexd.variables)
#print(bomexd.variables.values())
#########################################################################
# Define variables and importing files
bomexd = Dataset("bomex.default.0000000.nc","r")
area  = bomexd.variables["area"][:,:] # fractional area contained in mask
areah = bomexd.variables['areah'][:,:]
time_t=bomexd.variables['time'][:] 

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
vgrad=bomexd.variables['vgrad'][:,:]
u=bomexd.variables['u'][:,:] # vel in u direction
v=bomexd.variables['v'][:,:] # vel in v direction
w=bomexd.variables['w'][:,:] # vertical velocities



b=bomexd.variables['b'][:,:] # buoyancy
rho=bomexd.variables['rho'][:,:] # full level density
rhoh=bomexd.variables['rhoh'][:,:] # half level density
phydro=bomexd.variables['phydro'][:,:] # full level hydrostatic pressure


bomexql=Dataset("bomex.ql.0000000.nc","r")
qlfrac1=bomexql.variables['qlfrac'][:,:]
z1=bomexql.variables['z'][:]
zh1=bomexql.variables['zh'][:]


bomexqlcore=Dataset("bomex.qlcore.0000000.nc","r")
qlfrac2=bomexqlcore.variables['qlfrac'][:,:]
z2=bomexqlcore.variables['z'][:]
zh2=bomexqlcore.variables['zh'][:]


bomextrack=Dataset("l.0001800.track.nc","r") 
cv=bomextrack.variables['cv'][:] # cloud cover by volume 
cp=bomextrack.variables['cp'][:] # cloud cover by area
area_proj=bomextrack.variables['area_proj'][:] # projected cloud area
cloudsize=np.sqrt(area_proj)
area_proj_conn=bomextrack.variables['area_proj_conn'][:] # connected projected cloud area
ht=bomextrack.variables['ht'][:] # height
zt=bomextrack.variables['z'][:]
xt=bomextrack.variables['x'][:]
yt=bomextrack.variables['y'][:]
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

nz=zt.size;nx=xt.size;ny=yt.size; # sizes of dimensions




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



bomexw=Dataset('w.nc','r')

###########################################################################

def pyth(u,v):
    return np.sqrt(u*u+v*v)
    

def r1(x):
    return 1/(1+abs(x))

def r2(x):
    return 1/(1+x*x)
    
def rpow(x):
    return 1/(1+(abs(x))**(0.2))


def rb(u,w,h,b):
    return 1/(1+abs((u/w)*(h/b)))

def rbpow(u,w,h,b):
    return 1/(1+abs((u**0.2/w**0.2)*(h**0.2/b**0.2)))


def rexp(x):
    return np.exp(-abs(x))



def rb1(grad,ht,w,one,width):
    return 1/(abs(grad)*ht*(1/width*np.maximum(abs(w),one))+1)

def ra0(grad,ht,w,one):
    multiple=(abs(grad)*ht)/np.maximum(abs(w),one)
    return (1+multiple)**(-1)



def ra(grad,ht,w,ct,cb,one):
    multiple=(abs(grad)*ht)/np.maximum(abs(w),one)
    multiplicand=np.maximum(abs(ct-cb),one)/ct
    return (1+multiple*multiplicand)**(-1)




####################################
    
filenames=[bomexd]#, ricod, armd]

for file in filenames:
    zt=file.variables['z'][:]
    zh=file.variables['zh'][:]
    time_t=file.variables['time'][:]
    ugrad=file.variables['ugrad'][:,:]
    vgrad=file.variables['vgrad'][:,:]
    u=file.variables['u'][:,:]
    v=file.variables['v'][:,:]
    w=file.variables['w'][:,:]
    
    
    

    plt.figure()
    for t in range(time_t.size):
        unew=u[t,:]
        vnew=v[t,:]
        uv=pyth(unew,vnew)
        grad=np.gradient(uv)
        
        
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        #plt.figure()
        #plt.plot(zt, ugrad2(zt))
        #r=1/(1+abs(ugrad2(zt)))
        #plt.plot(ugrad2(zt),r1(ugrad2(zt)))
        
        #plt.plot(grad,r1(grad))
        #plt.title('r vs. gradient ')
        plt.plot(1/r1(grad),grad)
        plt.title('gradient vs. 1/r')
        
        
        ##
        gradint=interp1d(zt,grad,axis=0)
        ztnew=np.linspace(zt[0],zt[-1],int(overlap_ratio.size))
        #plt.plot(gradint(ztnew),overlap_ratio,'o')
        ##
        
        bins=15
        plt.hist2d(1/overlap_ratio, gradint(ztnew), bins=bins,cmin=0.5)
        #colorbar = plt.colorbar()
        #colorbar.set_label('counts in bin')
   


            # using rb1
    plt.figure()
    for t in range(time_t.size):
        unew=u[t,:]
        vnew=v[t,:]
        uv=pyth(unew,vnew)
        grad=np.gradient(uv)
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        gradint=interp1d(zt,grad,axis=0)
        wnew1=w[t,:]
        wnew2=interp1d(zh,wnew1,axis=0)
        #rb=1/(1+(ugrad2(ztnew)/1)*(ht/cb))
        
        ztnew=np.linspace(zt[0],zt[-1],int(ht.size))
        one=np.ones(ztnew.size) # let w=1 m/s, a way avoid dividing by zero
        
        plt.plot(1/rb1(gradint(ztnew),ht,wnew2(ztnew),one,np.sqrt(cv)),gradint(ztnew),'o')
        plt.title('gradient vs. 1/ratio')
        bins=15
        plt.hist2d(1/overlap_ratio, gradint(ztnew),bins=bins,cmin=0.5)
        #plt.plot(gradint(ztnew),rb1(gradint(ztnew),ht,wnew2(ztnew),one,np.sqrt(cv)),'o')
        #plt.title('ratio vs. grad')
        


                # using ra
    plt.figure()
    for t in range(time_t.size):
        unew=u[t,:]
        vnew=v[t,:]
        uv=pyth(unew,vnew)
        grad=np.gradient(uv)
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        gradint=interp1d(zt,grad,axis=0)
        wnew1=w[t,:]
        wnew2=interp1d(zh,wnew1,axis=0)
        #rb=1/(1+(ugrad2(ztnew)/1)*(ht/cb))
        
        ztnew=np.linspace(zt[0],zt[-1],int(ht.size))
        one=np.ones(ztnew.size) # let w=1 m/s, a way avoid dividing by zero
        
        plt.plot(1/ra(gradint(ztnew),ht,wnew2(ztnew),ct,cb,one),gradint(ztnew),'o')
        plt.title('gradient vs. 1/ratio')
        bins=15
        plt.hist2d(1/overlap_ratio, gradint(ztnew),bins=bins,cmin=0.5)
        #plt.xlim([1,1.8])
        ### this is for a better intuitive sense of what is going on
        #plt.plot(gradint(ztnew),ra0(gradint(ztnew),ht,wnew2(ztnew),one),'o')
        #plt.title('ratio vs. grad')
        #plt.hist2d(gradint(ztnew),overlap_ratio,bins=bins,cmin=0.5)
    
  ############################################################ 
# shear vs. 1/ overlap with hist2d grid
"""
    plt.figure()
    for t in range(time_t.size):
        unew=u[t,:]
        vnew=v[t,:]
        uv=pyth(unew,vnew)
        grad=np.gradient(uv)
        
        
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        #plt.figure()
        #plt.plot(zt, ugrad2(zt))
        #r=1/(1+abs(ugrad2(zt)))
        #plt.plot(ugrad2(zt),r1(ugrad2(zt)))
        
        #plt.plot(grad,r1(grad))
        #plt.title('r vs. gradient ')
        plt.plot(1/r1(grad),grad)
        plt.title('gradient vs. 1/r')
        
        
        ##
        gradint=interp1d(zt,grad,axis=0)
        ztnew=np.linspace(zt[0],zt[-1],int(overlap_ratio.size))
        #plt.plot(gradint(ztnew),overlap_ratio,'o')
        ##
        
        bins=15
        plt.hist2d(1/overlap_ratio, gradint(ztnew), bins=bins,cmin=0.5)
        #colorbar = plt.colorbar()
        #colorbar.set_label('counts in bin')
   


            # using rb(u,w,h,b)
    plt.figure()
    for t in range(time_t.size):
        unew=u[t,:]
        vnew=v[t,:]
        uv=pyth(unew,vnew)
        grad=np.gradient(uv)
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        gradint=interp1d(zt,grad,axis=0)
        wnew1=w[t,:]
        wnew2=interp1d(zh,wnew1,axis=0)
        one=np.ones(ztnew.size) # let w=1 m/s, a way avoid dividing by zero
        #rb=1/(1+(ugrad2(ztnew)/1)*(ht/cb))
        
        ztnew=np.linspace(zt[0],zt[-1],int(ht.size))
        plt.plot(1/ra(gradint(ztnew),ht,wnew2(ztnew),one,np.sqrt(cv)),gradint(ztnew),'o')
        plt.title('gradient vs. 1/ratio')
        bins=15
        plt.hist2d(1/overlap_ratio, gradint(ztnew),bins=bins,cmin=0.5)

"""

###################################################################################
     
"""
# seems like vert updraft vel has little affect on ratio
    plt.figure()
    for t in range(time_t.size):
        wnew1=w[t,:]
        
        
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        #plt.figure()
        #plt.plot(zt, ugrad2(zt))
        #r=1/(1+abs(ugrad2(zt)))
        #plt.plot(ugrad2(zt),r1(ugrad2(zt)))
        
        
        
        
        ##
        wint=interp1d(zh,wnew1,axis=0)
        ztnew=np.linspace(zt[0],zt[-1],int(overlap_ratio.size))
        #plt.plot(gradint(ztnew),overlap_ratio,'o')
        ##
        
        bins=15
        plt.hist2d(wint(ztnew), overlap_ratio, bins=bins,cmin=10.5)
        #plt.plot(wint(ztnew),overlap_ratio,'o')
        #colorbar = plt.colorbar()
        #colorbar.set_label('counts in bin')
"""


"""
        # using rb(u,w,h,b)
    plt.figure()
    for t in range(time_t.size):
        unew=u[t,:]
        vnew=v[t,:]
        uv=pyth(unew,vnew)
        grad=np.gradient(uv)
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        gradint=interp1d(zt,grad,axis=0)
        wnew1=w[t,:]
        wnew2=interp1d(zh,wnew1,axis=0)
        one=np.ones(ztnew.size) # let w=1 m/s, a way avoid dividing by zero
        #rb=1/(1+(ugrad2(ztnew)/1)*(ht/cb))
        
        ztnew=np.linspace(zt[0],zt[-1],int(ht.size))
        plt.plot(gradint(ztnew),rb(gradint(ztnew),one,ht,cb),'o')
        plt.title('rb vs. gradient ')
        bins=15
        plt.hist2d(gradint(ztnew), overlap_ratio, bins=bins,cmin=0.5)
          
        # plt.plot(zt,wnew2(zt)) # w is small ~0
"""    

    
"""  
       # plotting grad vs. w 
    plt.figure()
    for t in range(time_t.size):
        unew=u[t,:]
        vnew=v[t,:]
        uv=pyth(unew,vnew)
        grad=np.gradient(uv)
        #ugrad1=ugrad[t,:]
        #ugrad2=interp1d(zh,ugrad1)
        gradint=interp1d(zt,grad,axis=0)
        wnew1=w[t,:]
        wnew2=interp1d(zh,wnew1,axis=0)
        #one=np.ones(ztnew.size) # let w=1 m/s, a way avoid dividing by zero
        #rb=1/(1+(ugrad2(ztnew)/1)*(ht/cb))
        #ztnew=np.linspace(zt[0],zt[-1],int(ht.size))
        #plt.plot(gradint(ztnew),rb(gradint(ztnew),one,ht,cb),'o')
        #plt.title('rb vs. gradient ')
        plt.plot(wnew2(zt),gradint(zt))
        plt.title('grad vs. w')
"""
################################################################################
"""
    plt.figure()
    for t in range(time_t.size):
        vgrad1=vgrad[t,:]
        vgrad2=interp1d(zh,vgrad1)
        #plt.figure()
        #plt.plot(zt, ugrad2(zt))
        r=1/(1+abs(vgrad2(zt)))
        plt.plot(vgrad2(zt),r)
        plt.title('r vs. gradient of v')
"""




"""
    for t in range(time_t.size):
        plt.figure()
        ugrad1=ugrad[t,:]
        ugrad2=interp1d(zh,ugrad1)
        wnew1=w[t,:]
        wnew2=interp1d(zh,wnew1)
        plt.plot(wnew2(zt),ugrad2(zt))
        plt.title('gradient of u vs. w')
"""


"""
plt.figure()
plt.plot(cb,ht) # cb= width of cld base?
plt.title('ht vs. cb')

plt.figure()
plt.plot(cb,ct) 
plt.title('ct vs. cb')
plt.show()
"""

"""
# r(u,w,ht,cb)
ztnew=np.linspace(zt[0],zt[-1],int(ht.size))

plt.figure()
for t in range(time_t.size):
    ugrad1=ugrad[t,:]
    ugrad2=interp1d(zh,ugrad1)
    wnew1=w[t,:]
    wnew2=interp1d(zh,wnew1)
    
    plt.plot(wnew2(ztnew),ugrad2(ztnew))
    
    #rb=1/(1+(ugrad2(ztnew)/wnew2(ztnew))*(ht/cb))
    
    #plt.plot(ugrad2(ztnew),rb)
    #plt.title('rb vs. gradient of u')
"""

"""
#   r(u,w,ht,cb) with w=1
ztnew=np.linspace(zt[0],zt[-1],int(ht.size))

plt.figure()
for t in range(time_t.size):
    unew=u[t,:]
    vnew=v[t,:]
    uv=pyth(unew,vnew)
    grad=np.gradient(uv)
    #ugrad1=ugrad[t,:]
    #ugrad2=interp1d(zh,ugrad1)
    zh=bomexd.variables['zh'][:] 
    wnew1=w[t,:]
    wnew2=interp1d(zh,wnew1)
    
    #rb=1/(1+(ugrad2(ztnew)/1)*(ht/cb))
    
    plt.plot(grad,rb(grad,wnew2(z),ht,cb))
    plt.title('rb vs. gradient of u')

"""

"""
# plots of overlap_ratio (track file) vs. grad of u
for t in range(time_t.size):
    ugrad1=ugrad[t,:]
    ugrad2=interp1d(zh,ugrad1)
    ztnew=np.linspace(zt[0],zt[-1],int(overlap_ratio.size))
    plt.figure()
    #plt.plot(ugrad2(ztnew),overlap_ratio)
    bins=15
    plt.hist2d(ugrad2(ztnew), overlap_ratio, bins=bins,cmin=0.5)
    colorbar = plt.colorbar()
    colorbar.set_label('counts in bin')
    plt.title('ratio vs. gradient of u')

 """








#############################################################################
# shows that ct>=cb
# plt.plot(cb,ct,'o');plt.plot(cb,cb,'-');None
# ct -cb incr as r decr *(inv rel)
#plt.plot(ct-cb,overlap_ratio,'o');None
# linear relationship 
#plt.plot(ct-cb,ht,'o');None
# inv relationship
# plt.plot(ht,overlap_ratio,'o');None
# inv relationship
# plt.plot(area_proj,overlap_ratio,'o');None





####################################################################
#savfig('filename.pdf')
#savfig('filename.png')



end= time.time()
print('Run Time in Seconds:', end-start)













