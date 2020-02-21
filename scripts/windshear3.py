#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 12:00:15 2019

@author: anthonys
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 15:33:07 2019

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
"""
armd= Dataset("arm.default.0000000.nc","r")
armql=Dataset("arm.ql.0000000.nc","r")
armqlcore=Dataset("arm.qlcore.0000000.nc","r")
ricod= Dataset("rico.default.0000000.nc","r")
ricoql=Dataset("rico.ql.0000000.nc","r")
ricoqlcore=Dataset("rico.qlcore.0000000.nc","r")
"""
#data = pkgutil.get_data("conf", "fielddump.ql.00.05.track.nc")


#print(bomexd.variables.keys())
#print(bomexd.variables['area'])
#print(bomexd.variables['area'][:][:])
#print(bomexd.variables)
#print(bomexd.variables.values())
#########################################################################
# Define variables and importing files
"""
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

cv1=bomextrack1.variables['cv'][:]
cp1=bomextrack1.variables['cp'][:]


bomextrack2=Dataset("l.0005400.track.nc","r") 

cv2=bomextrack.variables['cv'][:]
cp2=bomextrack.variables['cp'][:]
area_proj2=bomextrack.variables['area_proj'][:]
cloudsize2=np.sqrt(area_proj)
ht2=bomextrack.variables['ht'][:]




bomexw=Dataset('w.nc','r')
bomexqt=Dataset('/data/bomex/qt.nc','r')

"""
################################################################################
# importing files
bomexd = Dataset("/data/bomex/bomex.default.0000000.nc","r")
bomexql = Dataset("/data/bomex/bomex.ql.0000000.nc","r")
bomexqlcore = Dataset("/data/bomex/bomex.qlcore.0000000.nc","r")
bomextrack18 = Dataset('/data/bomex/l.0001800.track.nc','r')
bomextrack36 = Dataset('/data/bomex/l.0003600.track.nc','r')
bomextrack54 = Dataset('/data/bomex/l.0005400.track.nc','r')
bomextrack72 = Dataset('/data/bomex/l.0007200.track.nc','r')



armd = Dataset("/data/arm/arm.default.0000000.nc","r")
armql = Dataset("/data/arm/arm.ql.0000000.nc","r")
armqlcore = Dataset("/data/arm/arm.qlcore.0000000.nc","r")
armtrack108 = Dataset('/data/arm/l.0010800.track.nc','r')
armtrack126 = Dataset('/data/arm/l.0012600.track.nc','r')
armtrack144 = Dataset('/data/arm/l.0014400.track.nc','r')
armtrack162 = Dataset('/data/arm/l.0016200.track.nc','r')


ricod = Dataset("/data/rico/rico.default.0000000.nc","r")
ricoql = Dataset("/data/rico/rico.ql.0000000.nc","r")
ricoqlcore = Dataset("/data/rico/rico.qlcore.0000000.nc","r")
ricotrack36 = Dataset('/data/rico/l.0003600.track.nc','r')
ricotrack72 = Dataset('/data/rico/l.0007200.track.nc','r')
ricotrack108 = Dataset('/data/rico/l.0010800.track.nc','r')
ricotrack144 = Dataset('/data/rico/l.0014400.track.nc','r')




###########################################################################

def pyth(u,v):
    return np.sqrt(u*u+v*v)
    

def r1(x):
    return 1/(1+abs(x))


def rb1(grad,ht,w,one,width):
    return 1/(abs(grad)*ht*(1/width*np.maximum(abs(w),one))+1)

def ra0(grad,ht,w,one):
    multiple=(abs(grad)*ht)/np.maximum(abs(w),one)
    return (1+multiple)**(-1)

def ra1(grad,ht,w,ct,cb,one):
    multiple=(abs(grad)*ht)/np.maximum(abs(w),one)
    multiplicand=1/ct
    return (1+multiple*multiplicand)**(-1)

def ra2(grad,ht,w,ct,cv,cb,one):
    multiple=(abs(grad)*ht)/np.maximum(abs(w),one)
    multiplicand=np.maximum(abs(ct-cb),one)/(cv)**(0.5)
    return (1+multiple*multiplicand)**(-1)

def ra(grad,ht,w,ct,cb,one):
    multiple=(abs(grad)*ht)/np.maximum(abs(w),one)
    multiplicand=np.maximum(abs(ct-cb),one)/ct
    return (1+multiple*multiplicand)**(-1)


####################################
    



filenames=[bomexd]#, ricod, armd]
bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144]




for file in filenames:
    zt=file.variables['z'][:]
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
    elif file == ricod:
        for file1 in ricofilenames:
            ht=file1.variables['ht'][:]
            cb=file1.variables['cb'][:]
            ct=file1.variables['ct'][:]
            cv=file1.variables['cv'][:]
            cp=file1.variables['cp'][:]
            overlap_ratio=file1.variables['chr'][:]
            area_proj=file1.variables['area_proj'][:]
    elif file == armd:
        for file1 in armfilenames:
            ht=file1.variables['ht'][:]
            cb=file1.variables['cb'][:]
            ct=file1.variables['ct'][:]
            cv=file1.variables['cv'][:]
            cp=file1.variables['cp'][:]
            overlap_ratio=file1.variables['chr'][:]
            area_proj=file1.variables['area_proj'][:]
            
        



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
                #plt.figure()
                #plt.plot(1/ra(gradint(ztnew),ht,wnew2(ztnew),ct,cb,one),gradint(ztnew),'-')
                #plt.title('gradient vs. 1/ratio')
                bins=15
                #plt.hist2d(1/overlap_ratio, gradint(ztnew),bins=bins,cmin=0.5)
                #plt.xlim([1,1.8])
                ### this is for a better intuitive sense of what is going on
                plt.plot(ra(gradint(ztnew),ht,wnew2(ztnew),ct,cb,one),gradint(ztnew),'-')
                plt.title('grad vs. ratio')
                ##plt.hist2d(gradint(ztnew),overlap_ratio,bins=bins,cmin=0.5)
                #plt.hist2d(ra(gradint(ztnew),ht,wnew2(ztnew),ct,cb,one),ct-cb,bins=bins,cmin=0.5)
                #plt.title('ct vs. ratio')
                
        
          
        
        """
              # using rb1
               #plt.figure()
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
                plt.figure()
                #plt.plot(1/rb1(gradint(ztnew),ht,wnew2(ztnew),one,np.sqrt(cv)),gradint(ztnew),'-')
                #plt.title('gradient vs. 1/ratio')
                bins=15
                #plt.hist2d(1/overlap_ratio, gradint(ztnew),bins=bins,cmin=0.5)
                plt.plot(rb1(gradint(ztnew),ht,wnew2(ztnew),one,np.sqrt(cv)),gradint(ztnew),'-')
                plt.title('grad vs. ratio')
        """
        
        
        
        
        
        """
            # using ra0
            #plt.figure()
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
                plt.figure()
                #plt.plot(1/ra0(gradint(ztnew),ht,wnew2(ztnew),one),gradint(ztnew),'-')
                #plt.title('gradient vs. 1/ratio')
                bins=15
                #plt.hist2d(1/overlap_ratio, gradint(ztnew),bins=bins,cmin=0.5)
                #plt.xlim([1,1.8])
                ### this is for a better intuitive sense of what is going on
                plt.plot(ra0(gradint(ztnew),ht,wnew2(ztnew),one),gradint(ztnew),'-')
                plt.title('grad vs. ratio')
        
        """



####################################################################
#savfig('filename.pdf')
#savfig('filename.png')



end= time.time()
print('Run Time in Seconds:', end-start)



