#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:47:58 2019

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


plt.rcParams.update({'font.size': 16})

start=time.time()


#############################################################################

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

####################################
    
# script

epsilon1=50    #arbitrarily set, most outliers are > 200, most nonoutliers clds are about 0

kt=0;
file_numb=-1
file1_numb=-1



### rico828

height=np.load('height_no_shear_rico828.npy')
volume=np.load('volume_no_shear_rico828.npy')
projar=np.load('projarea_no_shear_rico828.npy')
overlap=np.load('overlap_no_shear_rico828.npy')
areaz=np.load('area_z_ratio_rico828.npy')
areazg=np.load('rico_area_gmean82800.npy')
overlap_convex=np.load('overlap_convex_rico82800.npy')

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

conditional_height=0

file1_numb=-1

### accessing multiple datasets easily
#Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
#Bfilenames = Afilenames[5:]
Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
Bfilenames = Afilenames[22:23]
#Bfilenames = Afilenames[55:56]

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
            
    binsize=0.5
    condition_vec= np.arange(binsize,max(ht/l)+binsize,binsize)  #, 500, 100, 50] # to go through multiple conditions quickly
    h_shear = np.zeros(condition_vec.size)
    h_shape = np.zeros(condition_vec.size)
    h_turb = np.zeros(condition_vec.size)
    h_ACT = np.zeros(condition_vec.size)
    h_RES = np.zeros(condition_vec.size)
    m=-1
    
    
    for conditional_height in condition_vec:
        
        m=m+1
        #conditional_height=1500
        index=np.where( (conditional_height-binsize < ht/l)  & (ht/l <= conditional_height)  ) # indices of where condition holds true in ht vector
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
        
        #fractalf = np.minimum(((3*np.pi) / 4)*np.ones(overlap_ratio.size),1/overlap_ratio) 
        
        #fractalf = np.minimum( (3*np.pi /4)*np.ones(overlap_ratio.size), 1/overlap_ratio )
        #fractalf = np.minimum(fractal*np.ones(overlap_ratio.size),1/overlap_ratio)
        
        
        r_shear= 1/overlap_ratio[index[0]] -  1/overlap[index[0]]  # index variables with index_greater[0]
        r_shear[r_shear<0]=0
        r_area=1/areaz[index[0]] -1
        #r_turb= (np.ones(index_size)*fractal) -1
        r_turb = fractalf[index[0]]
        actual = 1/overlap_ratio[index[0]] 
        
        r_total =  ( r_shear + r_area + r_turb ) +1  # correction: need to add 1
        
        residual = actual - r_total 
        #residual[residual<0]=0
        
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
        """
        print('average r  is %.2f' % r3_avg, ', r std is %.2f' % r3_std)
        print('average r_shear  is %.2f' % S3_avg, ', r_shear std is %.2f' % S3_std)
        print('average r_area percentage is %.2f' % A3_avg, ', r_area std is %.2f' % A3_std)
        print('average r_turbulence percentage is %.2f' % T3_avg, ', r_turbulence std is %.2f' % T3_std)
        """
        
        
        ######################################################
        
      
    
            
            


            
h_shear=np.nan_to_num(h_shear)
h_shape=np.nan_to_num(h_shape)
h_turb=np.nan_to_num(h_turb)          
h_ACT=np.nan_to_num(h_ACT) 
h_RES=np.nan_to_num(h_RES)   

h_ACT=h_ACT -1
h_ACT[h_ACT<0]=0         
##########################################################################################



"""
h_shape=  ( np.load('h_shape_rico828.npy') + np.load('h_shape_rico2016.npy')  ) /2
h_shear=  ( np.load('h_shear_rico828.npy') + np.load('h_shear_rico2016.npy')  ) /2
h_turb=  ( np.load('h_turb_rico828.npy') + np.load('h_turb_rico2016.npy')  ) /2
"""

N= int(condition_vec.size)

ind = np.arange(N)    # the x locations for the groups
width = 0.8 #1/N     # the width of the bars: can also be len(x) sequence

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
plt.xlabel('Cloud Aspect Ratio')
plt.ylabel('Average Overlap')
plt.title('Average Contributions of Factors')
#plt.xticks(ind, '%' %[condition_vec])
plt.xticks([0 , 9 , 19 ], ('0.5', '5', '10' ) )
plt.yticks(np.arange(0, 10, 1))
#plt.legend((p1[0], p2[0],p3[0],p4[0]), ('Turbulence','Shear', 'Area Variability','Residual'), loc='lower right', bbox_to_anchor=(1.40, 0.5),fontsize='small',ncol=1)
plt.legend((p1[0], p2[0],p3[0]), ('Turbulence','Shear', 'Area Variability'), loc='lower right', bbox_to_anchor=(1.50, 0.5),fontsize='small',ncol=1)

plt.plot(h_ACT,'ko-',linewidth=2)

            
#plt.savefig('bomexAll_factor_average1_052219.eps', dpi=300, bbox_inches='tight')          
            
#ht_bin=ind*100 +100
ht_bin=ind*0.5 +0.5

###################################################################################
##################################################################################
#### curve fit for shear


"""
Z=np.where(h_shear == 0)
h_shear=np.delete(h_shear,Z)
ht_bin=np.delete(ht_bin,Z)


opt_s, cov_s = curve_fit(f_shear , ht_bin[2:], h_shear[2:])

c=opt_s
model_shear=f_shear(ht/l,c)


bins=1*dz;

plt.figure();
plt.hist2d(ht/l,model_shear,bins=[1*dz,1*dz],cmin=1, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
#plt.plot(horizontal,horizontal)
#plt.plot(vertical,vertical)


plt.title('Effect of Shear');
#plt.title('Effect of Shear');
#plt.title('Effect of Area Variability');
#plt.title('Effect of Shear, Area Var.');
#plt.title('Effect of Turbulence');
#plt.title('Relating Overlap and Height');

plt.xlabel('aspect ratio')
plt.ylabel('Cal Overlap')
#plt.xlabel('Cloud Overlap')
#plt.ylabel('Cloud Height')

#colorbar=plt.colorbar(extend='both');
colorbar=plt.colorbar()
colorbar.set_label('count in bin')
"""
##########################################################
##################################################################################
#### curve fit for area



Z=np.where(h_shape == 0)
h_shape=np.delete(h_shape,Z)
ht_bin=np.delete(ht_bin,Z)


opt_a, cov_a = curve_fit(f_area , ht_bin[2:], h_shape[2:])

c=opt_a
model_shape=f_area(ht/l,c)


bins=1*dz;

plt.figure();
plt.hist2d(ht/l,model_shape,bins=[1*dz,1*dz],cmin=1, cmap='viridis'  , norm=colors.LogNorm()); #norm=colors.LogNorm()
#plt.plot(horizontal,horizontal)
#plt.plot(vertical,vertical)


plt.title('Effect of Area');
#plt.title('Effect of Shear');
#plt.title('Effect of Area Variability');
#plt.title('Effect of Shear, Area Var.');
#plt.title('Effect of Turbulence');
#plt.title('Relating Overlap and Height');

plt.xlabel('aspect ratio')
plt.ylabel('Cal Overlap')
#plt.xlabel('Cloud Overlap')
#plt.ylabel('Cloud Height')

#colorbar=plt.colorbar(extend='both');
colorbar=plt.colorbar()
colorbar.set_label('count in bin')

##########################################################



######################################################
##replot stacked hist
plt.figure()
#ht_bin=ind*100 +100
ht_bin=ind*0.5 +0.5

h_shape=f_area(ht_bin,c)
#h_shear=f_shear(ht_bin,c)

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
plt.xlabel('Cloud Aspect Ratio')
plt.ylabel('Average Overlap')
plt.title('Average Contributions of Factors')
#plt.xticks(ind, '%' %[condition_vec])
plt.xticks([0 , 9 , 19 ], ('0.5', '5', '10' ) )
plt.yticks(np.arange(0, 10, 1))
#plt.legend((p1[0], p2[0],p3[0],p4[0]), ('Turbulence','Shear', 'Area Variability','Residual'), loc='lower right', bbox_to_anchor=(1.40, 0.5),fontsize='small',ncol=1)
plt.legend((p1[0], p2[0],p3[0]), ('Turbulence','Shear', 'Area Variability'), loc='lower right', bbox_to_anchor=(1.50, 0.5),fontsize='small',ncol=1)

plt.plot(h_ACT,'ko-',linewidth=2)

########################################################################################










###############################################
end= time.time()
print('Run Time in Seconds:', end-start)




