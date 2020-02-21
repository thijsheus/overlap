#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 18:22:28 2019

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
#import glob
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

plt.rcParams.update({'font.size': 16})

start=time.time()


#############################################################################

# definitions
def pyth(u,v):    # magnitude
    return np.sqrt(u*u+v*v)
    



################################################################################
# importing files
"""
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
"""

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

"""
filenames=[bomexd, ricod, armd]

bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72, bomextrack90, bomextrack108, bomextrack126, bomextrack144, bomextrack162, bomextrack180, bomextrack198, bomextrack216]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144, ricotrack180, ricotrack216, ricotrack252, ricotrack288, ricotrack324, ricotrack360, ricotrack396]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162, armtrack180, armtrack198, armtrack216, armtrack234, armtrack252, armtrack270, armtrack288]
"""
###########################################################################

####################################
    
# script
filenames=[ricod]
ricofilenames=[ricotrack2016]
#bomexfilenames=[bomextrack36]
armfilenames=[armtrack522]
#conditional_height=1000
epsilon1=50    #arbitrarily set, most outliers are > 200, most nonoutliers clds are about 0

kt=0;
file_numb=-1
file1_numb=-1

for file in filenames:
    file_numb = file_numb+1
    #zt=file.variables['z'][:]
    zh=file.variables['zh'][:]
    time_t=file.variables['time'][:]
    u=file.variables['u'][:,:]
    v=file.variables['v'][:,:]
    w=file.variables['w'][:,:]
    
    
            
    if file == ricod:
        for file1 in ricofilenames:
            file1_numb = file1_numb+1
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
            
            
            nrcloudarray = np.ma.getdata(nrcloud) # unmask array
            
            
            
            condition_vec=[1000]#, 500, 100, 50] # to go through multiple conditions quickly
            
            for conditional_height in condition_vec:
                
                #conditional_height=1500
                index_greater=np.where(ht >= conditional_height) # indices of where condition holds true in ht vector
                index_greater_size=index_greater[0].size
                ht_anvil=ht[index_greater[0]]   # taking the ht values according to indices above
                overlap_ratio_anvil=overlap_ratio[index_greater[0]];
                area_proj_anvil=area_proj[index_greater[0]]
                
              
                print('Clouds with height greater than',conditional_height,'meters')
                print('Clouds that satisfy the given condition: ',index_greater_size)
                
                ########################################
                
                
                
                ### rico828
                """
                height=np.load('height_no_shear_rico828.npy')
                volume=np.load('volume_no_shear_rico828.npy')
                projar=np.load('projarea_no_shear_rico828.npy')
                overlap=np.load('overlap_no_shear_rico828.npy')
                areaz=np.load('area_z_ratio_rico828.npy')
                """
                ### rico2016
                height=np.load('height_no_shear_rico2016.npy')
                volume=np.load('volume_no_shear_rico2016.npy')
                projar=np.load('projarea_no_shear_rico2016.npy')
                overlap=np.load('overlap_no_shear_rico2016.npy')
                areaz=np.load('area_z_ratio_rico2016.npy')
                model=np.load('model_rico2016.npy')
                
                ############################################################
                
                fractal=  (6*np.sqrt(3)) / (5*np.pi) 
                
                r_shear= 1/overlap_ratio[index_greater[0]] -  1/overlap[index_greater[0]]  # index variables with index_greater[0]
                r_shear[r_shear<0]=0
                r_area=1/areaz[index_greater[0]] -1
                r_turb= 1/(np.ones(index_greater_size)*fractal) -1
                
                r_total =  ( r_shear + r_area + r_turb ) +1  # correction: need to add 1
                
                
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
                
                S3_avg=np.mean(S3)#*100
                A3_avg=np.mean(A3)#*100
                T3_avg=np.mean(T3)#*100
                S3_std=np.std(S3)
                A3_std=np.std(A3)
                T3_std=np.std(T3)
                
                print('average r  is %.2f' % r3_avg, ', r std is %.2f' % r3_std)
                print('average r_shear  is %.2f' % S3_avg, ', r_shear std is %.2f' % S3_std)
                print('average r_area percentage is %.2f' % A3_avg, ', r_area std is %.2f' % A3_std)
                print('average r_turbulence percentage is %.2f' % T3_avg, ', r_turbulence std is %.2f' % T3_std)
                
                
                
                ######################################################
                
                print('----------------------------------------------------------------------------------')
                conditional_height=100
                index_less=np.where(ht <= conditional_height) # indices of where condition holds true in ht vector
                index_less_size=index_less[0].size
                ht_anvil=ht[index_less[0]]   # taking the ht values according to indices above
                overlap_ratio_anvil=overlap_ratio[index_less[0]];
                area_proj_anvil=area_proj[index_less[0]]
                print('Clouds with height less than',conditional_height,'meters')
                print('Clouds that satisfy the given condition: ',index_less_size)
                
                
                
                
                fractal=  (6*np.sqrt(3)) / (5*np.pi) 
                
                r_shear= 1/overlap_ratio[index_less[0]] - 1/overlap[index_less[0]]  # index variables with index_less[0]
                r_shear[r_shear<0]=0
                r_area=1/areaz[index_less[0]] -1
                r_turb= 1/(np.ones(index_less_size)*fractal) -1
                
                r_total =  ( r_shear + r_area + r_turb ) +1  
                
                
            
                S1=r_shear
                A1=r_area
                T1=r_turb
                R1=r_total
                r1_avg=np.mean(r_total)
                r1_std=np.std(r_total)
                
                S1_avg=np.mean(S1)#*100
                A1_avg=np.mean(A1)#*100
                T1_avg=np.mean(T1)#*100
                S1_std=np.std(S1)
                A1_std=np.std(A1)
                T1_std=np.std(T1)
                
                print('average r  is %.2f' % r1_avg, ', r std is %.2f' % r1_std)
                print('average r_shear  is %.2f' % S1_avg, ', r_shear std is %.2f' % S1_std)
                print('average r_area percentage is %.2f' % A1_avg, ', r_area std is %.2f' % A1_std)
                print('average r_turbulence percentage is %.2f' % T1_avg, ', r_turbulence std is %.2f' % T1_std)
            
                print('----------------------------------------------------------------------------------')
            
            
            
#############################################
                
print('----------------------------------------------------------------------------------')
#conditional_height=1500
index_middle=np.where( (100 < ht) & (ht < 500)) # indices of where condition holds true in ht vector
index_middle_size=index_middle[0].size
ht_anvil=ht[index_middle[0]]   # taking the ht values according to indices above
overlap_ratio_anvil=overlap_ratio[index_middle[0]];
area_proj_anvil=area_proj[index_middle[0]]
print('Clouds with height between 100 and 1000 meters')
print('Clouds that satisfy the given condition: ',index_middle_size)




fractal=  (6*np.sqrt(3)) / (5*np.pi) 

r_shear= 1/overlap_ratio[index_middle[0]] - 1/overlap[index_middle[0]]  # index variables with index_middle[0]
r_shear[r_shear<0]=0
r_area=1/areaz[index_middle[0]] -1
r_turb= 1/(np.ones(index_middle_size)*fractal) -1

r_total =  ( r_shear + r_area + r_turb ) +1  # correction: need to subtract 1 from each 



S2=r_shear
A2=r_area
T2=r_turb

r2_avg=np.mean(r_total)
r2_std=np.std(r_total)
                
S2_avg=np.mean(S2)#*100
A2_avg=np.mean(A2)#*100
T2_avg=np.mean(T2)#*100
S2_std=np.std(S2)
A2_std=np.std(A2)
T2_std=np.std(T2)

print('average r  is %.2f' % r2_avg, ', r std is %.2f' % r2_std)
print('average r_shear  is %.2f' % S2_avg, ', r_shear std is %.2f' % S2_std)
print('average r_area percentage is %.2f' % A2_avg, ', r_area std is %.2f' % A2_std)
print('average r_turbulence percentage is %.2f' % T2_avg, ', r_turbulence std is %.2f' % T2_std)

print('----------------------------------------------------------------------------------')


            
            
            
################################################################
                
N = 3
"""
h_shear = ( 7.67, 54.94 , 71.41)
h_shape = ( 5.17, 27.77 , 19.46)
h_turb = (87.19, 17.29,9.13)
"""

h_shear = ( S1_avg, S2_avg , S3_avg)
h_shape = ( A1_avg, A2_avg , A3_avg)
h_turb = (T1_avg, T2_avg, T3_avg)
ind = np.arange(N)    # the x locations for the groups
width = 1/N     # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind, h_shear, width, color='blue')
p2 = plt.bar(ind, h_shape, width, bottom=h_shear, color='purple' )
p3=plt.bar(ind,h_turb,width,bottom=(h_shape[0] + h_shear[0] , h_shape[1] + h_shear[1], h_shape[2]+h_shear[2]) , color='pink')

#plt.ylabel('Contribution Percentages')
#plt.title('Average Percentages')
plt.ylabel('Average Overlap')
plt.title('Average Contributions of Factors')
plt.xticks(ind, ( 'H < 100m' , '100m < H < 1000m', 'H > 1000m'))
plt.yticks(np.arange(0, 6, 1))
plt.legend((p1[0], p2[0],p3[0]), ('Shear', 'Area','Turbulence'), loc='lower right', bbox_to_anchor=(1.40, 0.5),fontsize='small',ncol=1)



#plt.savefig('rico2016_factor_average_042319.eps', dpi=300, bbox_inches='tight')        
            
            
            
            
            
            
            
###############################################
end= time.time()
print('Run Time in Seconds:', end-start)
