#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 13:21:07 2019

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
bomexd = Dataset("bomex.default.0000000.nc","r")
armd= Dataset("arm.default.0000000.nc","r")
ricod= Dataset("rico.default.0000000.nc","r")


#data = pkgutil.get_data("conf", "fielddump.ql.00.05.track.nc")


#print(bomexd.variables.keys())
#print(bomexd.variables['area'])
#print(bomexd.variables['area'][:][:])
#print(bomexd.variables)
#print(bomexd.variables.values())
#########################################################################
# Define variables and importing files

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

bomdaleprof=Dataset("profiles.001.nc","r")

bomexfdt4=Dataset("fielddump.ql.00.04.track.nc","r") # use these more, cvx
area2d=bomexfdt4.variables['area2d'][:,:]
ht=bomexfdt4.variables['ht'][:]
area_proj=bomexfdt4.variables['area_proj'][:]
areaproj_cvx=bomexfdt4.variables['areaproj_cvx'][:]
cb=bomexfdt4.variables['cb'][:]
ct=bomexfdt4.variables['ct'][:]
cv=bomexfdt4.variables['cv'][:]
cp=bomexfdt4.variables['cp'][:]
cfrac=bomexfdt4.variables['cfrac'][:]
overlap_ratio=bomexfdt4.variables['chr'][:]
cld_perimeters=bomexfdt4.variables['cld_perimeters'][:]

bomexfdu=Dataset("fielddump.u.001.nc")
xm=bomexfdu.variables['xm'][:]
yt=bomexfdu.variables['yt'][:]
zt=bomexfdu.variables['zt'][:]
time_t=bomexfdu.variables['time'][:]
u=bomexfdu.variables['u'][:,:,:,:]

bomexfdv=Dataset("fielddump.v.001.nc")
xt=bomexfdv.variables['xt'][:]
ym=bomexfdv.variables['ym'][:]
#zt=bomexfdv.variables['zt'][:]
#time=bomexfdv.variables['time'][:]
v=bomexfdv.variables['v'][:,:,:,:]

bomexfdw=Dataset("fielddump.w.001.nc")



bomexdalesamp=Dataset('sampling.001.nc','r')

bomexw=Dataset('w.nc','r')

###########################################################################


"""
ugrad1=np.zeros(int(zh.size))
vgrad1=np.zeros(int(zh.size))
#ugrad2=np.zeros(int(ht.size))
for i in range(zh.size):
    ugrad1[i]=np.sum(ugrad[:,i])/(int(time_t.size))
    vgrad1[i]=np.sum(vgrad[:,i])/(int(time_t.size))
    
#domain=np.linspace(zh[0],zh[-1],int(ht.size))

ugrad2=interp1d(zh,ugrad1)
vgrad2=interp1d(zh,vgrad1)


bins=15
plt.figure()
plt.hist2d(ht, ugrad2(ht), bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('ugrad vs. ht')

plt.figure()
plt.hist2d(ht, vgrad2(ht), bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('vgrad vs. ht')


plt.show()
"""

"""
# for default and l.track files
unew1=np.zeros(int(z.size))
vnew1=np.zeros(int(z.size))

for i in range(z.size):
    unew1[i]=np.sum(u[:,i])/(int(time_t.size))
    vnew1[i]=np.sum(v[:,i])/(int(time_t.size))
    


unew2=interp1d(z,unew1)
vnew2=interp1d(z,vnew1)


bins=15
plt.figure()
plt.hist2d(ht, unew2(ht), bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('u vs. ht')

plt.figure()
plt.hist2d(ht, vnew2(ht), bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('v vs. ht')


plt.show()
"""
######################################################################
# for profile and fielddump u,v files
unew1=np.sum(u,axis=0)/int(time_t.size)
unew2=np.sum(unew1,axis=1)/int(yt.size)
unew3=np.sum(unew2,axis=1)/int(xm.size)
# west-east vel
unew4=interp1d(zt,unew3)#,fill_value="extrapolate")


bins=15
plt.figure()
plt.hist2d(ht, unew4(ht), bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('u vs. ht')
    
"""
plt.figure()
plt.hist2d(area_proj, unew4(area_proj), bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('u vs. area_proj')
"""


# south-north vel
vnew1=np.sum(v,axis=0)/int(time_t.size)
vnew2=np.sum(vnew1,axis=1)/int(ym.size)
vnew3=np.sum(vnew2,axis=1)/int(xt.size)

vnew4=interp1d(zt,vnew3)#,fill_value="extrapolate")


#bins=15
plt.figure()
plt.hist2d(ht, vnew4(ht), bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('v vs. ht')



    



plt.figure()
plt.hist2d(zt, unew3, bins=bins,cmin=0.5)
colorbar = plt.colorbar()
colorbar.set_label('counts in bin')
plt.title('u vs. zt')

error=10
locu=np.array([0,0])
for i in range(1,nr.size):
    loccld=np.array([cb[i],cb[i]+ht[i]])
    for j in range(unew3.size):
        if abs(cb[i]-zt[j])<error and abs(cb[i]-zt[j])<abs(cb[i]-zt[j-1]):
            locu[0]=unew3[j]
        if abs(cb[i]+ht[i]-zt[j])<error and abs(cb[i]+ht[i]-zt[j])<abs(cb[i]+ht[i]-zt[j-1]):
            locu[1]=unew3[j]
    #locu=np.array([unew3[int(loccld[0])],unew3[int(loccld[1])]])
    plt.plot(loccld,locu,'-o')












plt.show()
    

end= time.time()
print('Run Time in seconds:', end-start)



























