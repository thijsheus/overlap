#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:40:13 2019

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



import plotly.plotly as py
import plotly.graph_objs as go


plt.rcParams.update({'font.size': 16})

start=time.time()


#############################################################################



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
#ricotrack144 = Dataset('/data/rico/l.0014400.track.nc','r')
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

lassotrack306 = Dataset('/data/lasso/sims/20160611_micro/l.0030600.track.nc','r')

"""
filenames=[bomexd, ricod, armd]

bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72, bomextrack90, bomextrack108, bomextrack126, bomextrack144, bomextrack162, bomextrack180, bomextrack198, bomextrack216]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144, ricotrack180, ricotrack216, ricotrack252, ricotrack288, ricotrack324, ricotrack360, ricotrack396]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162, armtrack180, armtrack198, armtrack216, armtrack234, armtrack252, armtrack270, armtrack288]
"""
###########################################################################
# definitions
def pyth(u,v):    # magnitude
    return np.sqrt(u*u+v*v)
    

#### definitions
def f_turb(h,parm):
    return (3*np.pi /4)* (h/(parm+h))

def f_shear(aspect,parm):
    return parm*aspect**1

def f_area(aspect,parm):
    return parm*aspect**1
#################################################################################
    
# script
filenames=[ricod]
ricofilenames=[ricotrack828]
bomexfilenames=[bomextrack360]
armfilenames=[armtrack522]
lassofilenames=[lassotrack306]
#conditional_height=1000
epsilon1=50    #arbitrarily set, most outliers are > 200, most nonoutliers clds are about 0

kt=0;
file_numb=-1
file1_numb=-1


"""
### rico828

height=np.load('height_no_shear_rico828.npy')
volume=np.load('volume_no_shear_rico828.npy')
projar=np.load('projarea_no_shear_rico828.npy')
overlap=np.load('overlap_no_shear_rico828.npy')
areaz=np.load('area_z_ratio_rico828.npy')
areazg=np.load('rico_area_gmean82800.npy')
overlap_convex=np.load('overlap_convex_rico82800.npy')
#shearh_vel = np.load('rico82800_shearh_vel.npy')
#sheara_vel = np.load('rico82800_sheara_vel.npy')
#wavg = np.load('rico82800_wavg.npy')
#shearTB_vel = np.load('rico82800_shearTB_vel.npy')
wavg1 = np.load('wavg1_rico82800.npy')
shearTB = np.load('shearTB_rico82800.npy')
shear_sum = np.load('shear_sum_rico82800.npy')
wz_max = np.load('wz_max_rico82800.npy')
wz_cb = np.load('wz_cb_rico82800.npy')
shift_distance = np.load('shift_distancerico82800.npy')
"""
"""
### rico2016
height=np.load('height_no_shear_rico2016.npy')
volume=np.load('volume_no_shear_rico2016.npy')
projar=np.load('projarea_no_shear_rico2016.npy')
overlap=np.load('overlap_no_shear_rico2016.npy')
areaz=np.load('area_z_ratio_rico2016.npy')
model=np.load('model_rico2016.npy')
model_inv=np.load('model_inv_rico2016.npy') 
model_opt=np.load('model_opt_rico2016.npy')
model_opt2=np.load('model_opt2_rico2016.npy')
overlap_convex=np.load('overlap_convex_rico201600.npy')
wavg1 = np.load('wavg1_rico201600.npy')
shearTB = np.load('shearTB_rico201600.npy')
shear_sum = np.load('shear_sum_rico201600.npy')
shift_distance = np.load('shift_distancerico201600.npy')
"""
"""
### bomex 342
height342=np.load('height_no_shear_bomex34200.npy')
volume342=np.load('volume_no_shear_bomex34200.npy')
projar342=np.load('projarea_no_shear_bomex34200.npy')
overlap342=np.load('overlap_no_shear_bomex34200.npy')
areaz342=np.load('area_z_ratio_bomex34200.npy')

### bomex 360
height=np.load('height_no_shear_bomex36000.npy')
volume=np.load('volume_no_shear_bomex36000.npy')
projar=np.load('projarea_no_shear_bomex36000.npy')
overlap=np.load('overlap_no_shear_bomex36000.npy')
areaz=np.load('area_z_ratio_bomex36000.npy')

### lasso 306
height=np.load('height_no_shear_lasso306.npy')
volume=np.load('volume_no_shear_lasso306.npy')
projar=np.load('projarea_no_shear_lasso306.npy')
overlap=np.load('overlap_no_shear_lasso306.npy')
areaz=np.load('area_z_ratio_lasso306.npy')
"""

### arm28800

height=np.load('height_no_shear_arm28800.npy')
volume=np.load('volume_no_shear_arm28800.npy')
projar=np.load('projarea_no_shear_arm28800.npy')
overlap=np.load('overlap_no_shear_arm28800.npy')
areaz=np.load('area_z_ratio_arm28800.npy')
overlap_convex=np.load('overlap_convex_arm28800.npy')
#shearh_vel = np.load('rico82800_shearh_vel.npy')
#sheara_vel = np.load('rico82800_sheara_vel.npy')
#wavg = np.load('rico82800_wavg.npy')
#shearTB_vel = np.load('rico82800_shearTB_vel.npy')
wavg1 = np.load('wavg1_arm28800.npy')
shearTB = np.load('shearTB_vel_arm28800.npy')
shear_sum = np.load('shear_sum_arm28800.npy')
shift_distance = np.load('shift_distancearm28800.npy')

"""
### arm41400

height=np.load('height_no_shear_arm41400.npy')
volume=np.load('volume_no_shear_arm41400.npy')
projar=np.load('projarea_no_shear_arm41400.npy')
overlap=np.load('overlap_no_shear_arm41400.npy')
areaz=np.load('area_z_ratio_arm41400.npy')
overlap_convex=np.load('overlap_convex_arm41400.npy')
#shearh_vel = np.load('rico82800_shearh_vel.npy')
#sheara_vel = np.load('rico82800_sheara_vel.npy')
#wavg = np.load('rico82800_wavg.npy')
#shearTB_vel = np.load('rico82800_shearTB_vel.npy')
wavg1 = np.load('wavg1_arm41400.npy')
shearTB = np.load('shearTB_vel_arm41400.npy')
shear_sum = np.load('shear_sum_arm41400.npy')
shift_distance = np.load('shift_distancearm41400.npy')
"""
### Bomex_All is Below

conditional_height=0

file1_numb=-1

### accessing multiple datasets easily
#Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
#Bfilenames = Afilenames[5:]
#Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
#Bfilenames = Afilenames[22:23]
#Bfilenames = Afilenames[55:56]
Afilenames = sorted(glob.glob('/data/arm/*.track.nc'))
#Bfilenames = Afilenames[1:]
Bfilenames = Afilenames[10:11]
#Bfilenames = Afilenames[17:18]

for file1 in Bfilenames:
    print(file1)

    data = Dataset(file1,'r')

    file1_numb = file1_numb+1
    
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
            
    """
    ### Bomex_All
    overlap_ratio=np.load('Bomex_overlap_ratio.npy')            
    height=np.load('Bomex_ht.npy')
    overlap=np.load('Bomex_shear.npy')
    areaz=np.load('Bomex_areaz.npy')
    ht=height
    cv=np.load('Bomex_cv.npy')
    cp=np.load('Bomex_cp.npy')
    overlap_convex=np.load('Bomex_hull.npy')
    #shearh_vel = np.load('Bomex_shearh_vel.npy')
    #sheara_vel = np.load('Bomex_sheara_vel.npy')
    #wavg = np.load('Bomex_wavg.npy')
    shearTB = np.load('Bomex_shearTB.npy')
    shear_sum = np.load('Bomex_shear_sum.npy')
    wz_max = np.load('Bomex_wz_max.npy')
    wz_cb = np.load('Bomex_wz_cb.npy')
    shift_distance = np.load('Bomex_shift_distance.npy')
    print('Bomex')
    """
    
    
    
    binsize=100
    condition_vec= np.arange(binsize,max(ht)+binsize,binsize)  #, 500, 100, 50] # to go through multiple conditions quickly
    h_shear = np.zeros(condition_vec.size)
    h_shape = np.zeros(condition_vec.size)
    h_turb = np.zeros(condition_vec.size)
    h_ACT = np.zeros(condition_vec.size)
    h_RES = np.zeros(condition_vec.size)
    h_total = np.zeros(condition_vec.size)
    m=-1
    

    ##### loop over bin hts
    for conditional_height in condition_vec:
        
        m=m+1
        #conditional_height=1500
        index=np.where( (conditional_height-binsize < ht)  & (ht <= conditional_height)  ) # indices of where condition holds true in ht vector
        index_size=index[0].size
        ht_sel=ht[index[0]]   # taking the ht values according to indices above
        overlap_ratio_sel=overlap_ratio[index[0]];
        #area_proj_sel=area_proj[index[0]]
        
        
        """
        print('Clouds with height greater than',conditional_height,'meters')
        print('Clouds that satisfy the given condition: ',index_size)
        """
        ########################################
        
        #areaz=areazg
        
                        
        ############################################################
        ###### calculating overlap due to contributors
        l=cv**0.5
        fract= (5*np.pi) / (6*np.sqrt(3))  -1
        fractal=  (5*np.pi) / (6*np.sqrt(3)) # -1
        fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio) -1
        
        shape=  (1/areaz) - 1
        shear=   (1/overlap_ratio) - (1/overlap)  
        #shear = abs(shear)
        shear[shear<0]=0
        
        overlap_convex=overlap_convex[overlap_convex>0]
        turb= (1/overlap_ratio[ht>100]) - (1/overlap_convex)
        turb[turb<0]=0
        fractalf[ht>100]=turb
        
        
        ## upper bound to turb
        lam =   (3*np.pi)/4  #np.pi * np.sqrt(3) / 2     #(15/8)*np.sqrt(3/8)*np.pi
        #C1=np.argwhere((1/overlap_ratio)>1.5*lam)
        keep=np.argwhere(fractalf > lam)
        fractalf[np.array(keep)] = lam
        
      
        
        FF = f_turb(ht,200);FF[ht<=100]=0;
        shear_sum = abs(shear_sum)
        shift_d = np.zeros(ht.size)
        shift_d[ht>50] = shift_distance
        
        shift_d = shift_d + shape*cp**0.5/2  ### draw the picture
        updraft = (shear_sum*ht)/shift_d
        FS = (shear_sum*ht)/(updraft*cp**0.5)
        FS[updraft<0.1] = 0
        #FS[updraft<0.01] = 0 
        
        #shear = FS ; fractalf = FF
        
        #r_shear= 1/overlap_ratio[index[0]] -  1/overlap[index[0]]  # index variables with index_greater[0]
        #r_shear[r_shear<0]=0
        r_shear=shear[index[0]]
        r_area=1/areaz[index[0]] -1
        #r_turb= (np.ones(index_size)*fractal) -1
        r_turb = fractalf[index[0]]
        actual = 1/overlap_ratio[index[0]] 
        
        #r_total =  ( r_shear + r_area + r_turb ) +1  
        
        ########
        
        #### empirical model
        #f_t=f_turb(ht,200); f_t[ht<=100]=0
        #f_s=f_shear(ht/l,0.43);f_s[ht/l<=1]=0
        #f_a=f_area(ht/l,0.2);f_a[ht/l<=1]=0
        
        #r_shear=f_s[index[0]]
        #r_area=f_a[index[0]]
        #r_turb=f_t[index[0]]
        
        
        r_total =  ( r_shear + r_area + r_turb ) +1  
        
        
        ########
        
        residual = actual - r_total 
        #residual[residual<0]=0
        ####### binning data
        """
        S=((r_shear-1)/r_total)
        A=((r_area-1)/r_total)
        T=((r_turb-1)/r_total)
        """
        S3=r_shear
        A3=r_area
        T3=r_turb
        
        r3_avg=np.mean(r_total)
        r3_std=np.std(r_total)
        ACT_avg=np.mean(actual)
        ACT_std=np.std(actual)
        RES_avg=np.mean(residual)
        RES_std=np.std(residual)
        
        S3_avg=np.mean(S3)#*100
        A3_avg=np.mean(A3)#*100
        T3_avg=np.mean(T3)#*100
        S3_std=np.std(S3)
        A3_std=np.std(A3)
        T3_std=np.std(T3)
        
        h_shear[m]=S3_avg
        h_shape[m]=A3_avg
        h_turb[m]=T3_avg
        h_ACT[m]=ACT_avg
        h_RES[m]=RES_avg
        h_total[m] = r3_avg 
        """
        print('average r  is %.2f' % r3_avg, ', r std is %.2f' % r3_std)
        print('average r_shear  is %.2f' % S3_avg, ', r_shear std is %.2f' % S3_std)
        print('average r_area percentage is %.2f' % A3_avg, ', r_area std is %.2f' % A3_std)
        print('average r_turbulence percentage is %.2f' % T3_avg, ', r_turbulence std is %.2f' % T3_std)
        """
        
        
        ######################################################
        
      
    
            
            


"""
h_shear=np.nan_to_num(h_shear)
h_shape=np.nan_to_num(h_shape)
h_turb=np.nan_to_num(h_turb)          
h_ACT=np.nan_to_num(h_ACT) 
h_RES=np.nan_to_num(h_RES)   
"""
h_ACT=h_ACT -1
h_ACT[h_ACT<0]=0         
##########################################################################################
###### stacked bar graph
#h_turb=f_turb(ht_bin,opt_t)
## change variable
h_total = h_shape + h_shear +h_turb 


N= int(condition_vec.size)

#N=np.count_nonzero(~np.isnan(h_shear))

ind = np.arange(N)    # the x locations for the groups
width = 0.8 #1/N     # the width of the bars: can also be len(x) sequence
width = 90
ind=ind*100 +100
arg_nan=np.argwhere(np.isnan(h_shape))
#ind=ind[~np.isnan(h_shear)]
#N=ind.size
ind = ind.astype('float')
ind[arg_nan] = np.nan 
"""
ind=ind[~np.isnan(ind)]
h_shear=h_shear[~np.isnan(h_shear)]
h_shape=h_shape[~np.isnan(h_shape)]
h_turb=h_turb[~np.isnan(h_turb)]
h_ACT = h_ACT[~np.isnan(h_ACT)]
"""
bottom2=np.zeros(ind.size)
bottom3=np.zeros(ind.size)
for i in range(ind.size):
    #i = int(i/100 - 1)
    bottom2[i]=h_turb[i]+h_shear[i]
    bottom3[i]=h_turb[i]+h_shear[i]+h_shape[i]

#ind=ind[~np.isnan(ind)]
#h_shear=h_shear[~np.isnan(h_shear)]
#h_shape=h_shape[~np.isnan(h_shape)]
#h_turb=h_turb[~np.isnan(h_turb)]

"""
mystring = []
for digit in ind:
    mystring.append(digit)
ind = np.arange(100,100*ind.size+100,100)
"""
p1= plt.bar(ind,h_turb,width , color='pink')
p2 = plt.bar(ind, h_shear, width, bottom=h_turb  ,color='blue')
p3 = plt.bar(ind, h_shape, width, bottom=bottom2, color='purple' )
#p4 = plt.bar(ind, h_RES, width, bottom=bottom3, color='orange' )

#plt.ylabel('Contribution Percentages')
#plt.title('Average Percentages')
plt.xlabel('Cloud Height')
plt.ylabel('Average Overlap')
plt.title('Average Contributions of Factors')
#plt.xticks(ind,mystring)
#plt.xticks(ind[0:1],condition_vec[0:1])
#plt.xticks([0 , 9 , 19 ], ('100', '1000', '2000' ) )
plt.yticks(np.arange(0, 20, 1))
#plt.legend((p1[0], p2[0],p3[0],p4[0]), ('Turbulence','Shear', 'Area Variability','Residual'), loc='lower right', bbox_to_anchor=(1.40, 0.5),fontsize='small',ncol=1)


pa=plt.plot(ind,h_ACT,'ko',linewidth=2)

plt.legend((p1[0], p2[0],p3[0],pa[0]), ('Turbulence','Shear', 'Area Variability','Actual Overlap'), loc='lower right', bbox_to_anchor=(1.50, 0.5),fontsize='small',ncol=1)

#plt.savefig('arm28800_factor_average1_emp1_061719.eps', dpi=300, bbox_inches='tight')          
            

############################################################
#### makes table of fraction of factor, put it in matrix form 
print('Area', '\t', 'Shear',  '\t' ,'Turbulence')
h_1=h_shape/h_total
h_2=h_shear/h_total
h_3=h_turb/h_total

for j in range(h_ACT.size):
    print( '%.3f' %h_1[j],'\t', '%.3f '%h_2[j], '\t', '%.3f'  %h_3[j] )
    
print('average percentages (%.3f , %.3f , %.3f)' %(np.nanmean(h_1), np.nanmean(h_2), np.nanmean(h_3)))

h_matrix = np.array([[h_1],[h_2],[h_3]])
h_matrix2 = np.zeros((h_1.size,4))
h_matrix2[:,0] = np.round(ind,3)
h_matrix2[:,1] = np.round(h_1,3)
h_matrix2[:,2] = np.round(h_2,3)
h_matrix2[:,3] = np.round(h_3,3)
###################################################################################
##################################################################################
#### curve fit for turbulence

ht_bin=ind*100 +100
"""
Z=np.where(h_turb == 0)
h_turb=np.delete(h_turb,Z)
ht_bin=np.delete(ht_bin,Z)


opt_t, cov_t = curve_fit(f_turb , ht_bin[1:], h_turb[1:])

c=opt_t
model_turb=f_turb(ht,c)


bins=1*dz;

plt.figure();
plt.hist2d(ht,model_turb,bins=[1*dz,1*dz],cmin=1, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
#plt.plot(horizontal,horizontal)
#plt.plot(vertical,vertical)


plt.title('Effect of Turbulence');
#plt.title('Effect of Shear');
#plt.title('Effect of Area Variability');
#plt.title('Effect of Shear, Area Var.');
#plt.title('Effect of Turbulence');
#plt.title('Relating Overlap and Height');

plt.xlabel('ht')
plt.ylabel('Cal Overlap')
#plt.xlabel('Cloud Overlap')
#plt.ylabel('Cloud Height')

#colorbar=plt.colorbar(extend='both');
colorbar=plt.colorbar()
colorbar.set_label('count in bin')

##########################################################
"""



"""
######################################################
##replot stacked hist
plt.figure()
ht_bin=ind*100 +100
h_turb=f_turb(ht_bin,c)


bottom2=np.zeros(N)
bottom3=np.zeros(N)
for i in range(N):
    bottom2[i]=h_turb[i]+h_shear[i]
    bottom3[i]=h_turb[i]+h_shear[i]+h_shape[i]

p1= plt.bar(ind,h_turb,width , color='pink')
p2 = plt.bar(ind, h_shear, width, bottom=h_turb  ,color='blue')
p3 = plt.bar(ind, h_shape, width, bottom=bottom2, color='purple' )
#p4 = plt.bar(ind, h_RES, width, bottom=bottom3, color='orange' )

#plt.ylabel('Contribution Percentages')
#plt.title('Average Percentages')
plt.xlabel('Cloud Height')
plt.ylabel('Average Overlap')
plt.title('Average Contributions of Factors')
#plt.xticks(ind, '%' %[condition_vec])
plt.xticks([0 , 9 , 19 ], ('100', '1000', '2000' ) )
plt.yticks(np.arange(0, 8, 1))
#plt.legend((p1[0], p2[0],p3[0],p4[0]), ('Turbulence','Shear', 'Area Variability','Residual'), loc='lower right', bbox_to_anchor=(1.40, 0.5),fontsize='small',ncol=1)
plt.legend((p1[0], p2[0],p3[0]), ('Turbulence','Shear', 'Area Variability'), loc='lower right', bbox_to_anchor=(1.50, 0.5),fontsize='small',ncol=1)

plt.plot(h_ACT,'ko-',linewidth=2)

########################################################################################
"""


###############################################
end= time.time()
print('Run Time in Seconds:', end-start)






















