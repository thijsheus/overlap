#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 11:51:46 2019

@author: anthonys
"""
"""
calculates overlap ratio base on ht of cld from chopping the cld at the top


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

from scipy.stats import spearmanr


start=time.time()


#############################################################################

# definitions
def pyth(u,v):    # magnitude
    return np.sqrt(u*u+v*v)
    



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
ricotrack144 = Dataset('/data/rico/l.0014400.track.nc','r')
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


filenames=[bomexd, ricod, armd]

bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72, bomextrack90, bomextrack108, bomextrack126, bomextrack144, bomextrack162, bomextrack180, bomextrack198, bomextrack216]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144, ricotrack180, ricotrack216, ricotrack252, ricotrack288, ricotrack324, ricotrack360, ricotrack396]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162, armtrack180, armtrack198, armtrack216, armtrack234, armtrack252, armtrack270, armtrack288]

###########################################################################

####################################
    
##### choose case and track file
# script
filenames=[ricod]
ricofilenames=[ricotrack828]
bomexfilenames=[bomextrack36]
armfilenames=[armtrack522]
conditional_height=1500

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
            
            
            
            ### take cloud that satisfy conditional ht
            
            #conditional_height=1500
            index_anvil=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
            index_anvil_size=index_anvil[0].size
            ht_anvil=ht[index_anvil[0]]   # taking the ht values according to indices above
            overlap_ratio_anvil=overlap_ratio[index_anvil[0]];
            area_proj_anvil=area_proj[index_anvil[0]]
            print('Clouds that satisfy the given condition: ',index_anvil_size)
            print('conditional height is',conditional_height,'meters')
            
            
            



############################################################################
            # calculate overlap of a cld that at some of its height chopped off at the top 
            
            
            index7=index_anvil[0] # to make easy on coder
            dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
            gridarea=dx*dy
            gridvol=dx*dy*dz
            nx=xt.size;ny=yt.size;nz=zt.size;
            
           
            cloud_numb=index7 +1 # index is off by 1 as array starts with zero
            #countnew1=count_elmts[index7] # dont need as vol is calculatd in for loop
            
            
            # percent of cloud ht should change countnew1=> change vol, ht, proj area
        
            #cloud_numb=np.array([533]) # choosing a particular cld
            #vert_res=25
            height_c=np.zeros(cloud_numb.size)
            overlap_c=np.zeros(cloud_numb.size)
            projected_area_c=np.zeros(cloud_numb.size)
            volume_c=np.zeros(cloud_numb.size)
            m=0;
            ###### choppong cld at a certain fractional ht
            #alpha=np.arange(1,0,-1) # just 1
            alpha=np.arange(1,0.35,-0.1) # fraction
            #alpha=np.array([0.7,0.5,0.3,0.1])
            
            #Dist=np.zeros((cloud_numb.size , alpha.size -1))
            #Dcurve=np.zeros((cloud_numb.size,alpha.size))  # multi dimensional array
            #Dcurve2=np.zeros((cloud_numb.size,alpha.size))  # multi dimensional array
            #Dcurve3=np.zeros((cloud_numb.size,alpha.size))  # multi dimensional array
            
            Total_dr=np.zeros((cloud_numb.size,alpha.size -1))  # multi dimensional array
            mu=np.zeros(cloud_numb.size)
            sigma=np.zeros(cloud_numb.size)
            Keeptrack_tdr=np.zeros((cloud_numb.size,3))  #make abitrarily big
            alpha_matrix=np.zeros((cloud_numb.size,alpha.size))
            overlap_changed_matrix=np.zeros((cloud_numb.size,alpha.size))
            
            ##### loop for clouds
            for e1 in cloud_numb:  # may take a while to run 
                
                
                location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1
                height_c[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
                #projected_area_c[m]=dx*dy*len(np.unique(location[2]+location[1]*nx)) #len of unique # of cells that get proj
                partial_heights=alpha*height_c[m] # fraction of cld height
                #height_diff=height_c[m]-partial_heights
                
                overlap_changed=np.zeros(partial_heights.size)
                projected_area_changed=np.zeros(partial_heights.size)
                vol_changed=np.zeros(partial_heights.size)
                part_loc= -1
                for partial in partial_heights:  ### loop over fractional hts of cld
                    
                    part_loc= part_loc +1 # to index later for diff ht levels
                    for i in range(location[0].size):
                        if dz*(location[0][i]-location[0][0]+1) <= partial: # finding all cells below chopped cld ht (take highest cell)
                            celltop=i
                        else:
                            break
                    #if partial - dz*(location[0][celltop]-location[0][0]+1) >0: 
                    difference = partial - dz*(location[0][celltop]-location[0][0]+1) # diff btw chopped cld ht and ht of clds w/ full cells
                    #print(difference == 0 )
                    #print(difference)
                    for j in range(celltop,location[0].size):  #### find top grid pts
                        if  location[0][celltop] +1  == location[0][j]: # to use partial cells
                            cellabove=j
                        elif location[0][celltop] +1  < location[0][j]:
                            break
                    #print('partial',partial)    
                    #print('top',celltop)
                    #print('above',cellabove)
                    #loc_loc0=-1
                    if partial == height_c[m]: # difference=0, so stop at celltop
                        proj_area_loop1=np.zeros(celltop +1)
                        vol_loop1=np.zeros(celltop +1)
                        
                        """
                        ##### 2d plot of cld
                        bins=dz
                        plt.figure()
                        plt.hist2d(dx*location[2][0:celltop+1],dz*location[0][0:celltop+1],bins=bins,cmin=0.5)
                        plt.xlabel('x');plt.ylabel('z')
                        plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                        colorbar = plt.colorbar()
                        colorbar.set_label('counts in bin')
                        plt.figure()
                        plt.hist2d(dy*location[1][0:celltop+1],dz*location[0][0:celltop+1],bins=bins,cmin=0.5)
                        plt.xlabel('y');plt.ylabel('z')
                        plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                        colorbar = plt.colorbar()
                        colorbar.set_label('counts in bin')
                        """
                        
                    else:
                        proj_area_loop1=np.zeros(cellabove +1)
                        vol_loop1=np.zeros(cellabove +1)
                        
                        """
                        ##### 2d plot of cld
                        bins=dz
                        plt.figure()
                        plt.hist2d(dx*location[2][0:cellabove+1],dz*location[0][0:cellabove+1],bins=bins,cmin=0.5)
                        plt.xlabel('x');plt.ylabel('z')
                        plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                        colorbar = plt.colorbar()
                        colorbar.set_label('counts in bin')
                        plt.figure()
                        plt.hist2d(dy*location[1][0:cellabove+1],dz*location[0][0:cellabove+1],bins=bins,cmin=0.5)
                        plt.xlabel('y');plt.ylabel('z')
                        plt.title( 'Cloud Number %s with height %.f m'%(cloud_numb[m],partial))
                        colorbar = plt.colorbar()
                        colorbar.set_label('counts in bin')
                        """
                                
                        
                    for loc0 in range(location[0].size):
                        
                        #loc_loc0 = loc_loc0 +1
                        if loc0  <= celltop : 
                            proj_area_loop1[loc0]=(location[2][loc0]+location[1][loc0]*nx)
                        elif celltop < loc0 and loc0 <=cellabove:
                            proj_area_loop1[loc0]=(location[2][loc0]+location[1][loc0]*nx)
                        if loc0  <= celltop: # full cell vol
                            vol_loop1[loc0]=(dx*dy*dz)   
                        elif celltop < loc0 and loc0 <=cellabove: # partial cell vol
                            vol_loop1[loc0]=(dx*dy*difference) 
                        else:
                            break
                        #if dz*loc0 < parital:
                    #proj_area_loop1=np.asarray(proj_area_loop1)
                    #print('loop1',proj_area_loop1.size)
                    projected_area_changed[part_loc]=dx*dy*len(np.unique(proj_area_loop1)) # counting cells w/ ! id #  & proj-ing
                    #vol_loop1=np.asarray(vol_loop1)  
                    vol_changed[part_loc]=sum(vol_loop1)
                    overlap_changed[part_loc]=vol_changed[part_loc]/(partial_heights[part_loc]*projected_area_changed[part_loc])
                    if partial == height_c[m]:
                        overlap_c[m]=overlap_changed[part_loc] # make sure calculated overlap is correct for full ht
                        volume_c[m]=vol_changed[part_loc]
                        projected_area_c[m]=projected_area_changed[part_loc]
                        
                    
                    ### dist btw adjacent pts
                    #if part_loc>0:
                     #   Dist[m,part_loc-1]= np.sqrt( (alpha[part_loc]-alpha[part_loc-1])**2 + (overlap_changed[part_loc]-overlap_changed[part_loc-1])**2  ) 
                     
    
                ###derivative at each pt
                #vol_grad=np.gradient(vol_changed)
                #projA_grad=np.gradient(projected_area_changed)
                #ht_grad=np.gradient(partial_heights)
                #Dcurve2[m]=(vol_grad/volume_c[m])-(ht_grad/height_c[m])-(projA_grad/projected_area_c[m])
                #overlap_changed_grad=np.gradient(overlap_changed)
                #alpha_grad=np.gradient(alpha)
                #Dcurve[m]=alpha_grad/overlap_changed_grad
                #Dcurve3[m]=np.gradient(alpha/overlap_changed)
                
                ## calculates change in r of a cld
                dr=np.zeros(alpha.size -1) 
                for i in range(alpha.size -1 ) :
                    dr[i]=  overlap_changed[i+1] - overlap_changed[i]    # having abs lowers the high pt of the interval
                    Total_dr[m,i]=dr[i]
                    
                
                
                alpha_matrix[m,:]=alpha
                overlap_changed_matrix[m,:]=overlap_changed
                Keeptrack_tdr[kt,:]=[file_numb,file1_numb,e1]
                kt=kt+1
                """
                mu[m]=np.mean(dr)
                sigma[m]=np.std(dr)
                
                ### prints index of matrix where change in r is outside the confidence interval
                for i in range(alpha.size -1 ):
                    if abs(Total_dr[m,i]) >= mu[m]+1.96*sigma[m]:
                        print('index',(m,i))  #remeber indices start at zero
                """
                ###putting alpha = partial_heights / height_c[m] b/c partial_heights = alpha * height_c[m]
                plt.plot(overlap_changed,alpha,'o-',label=str(cloud_numb[m])) 
                plt.title('normalized height vs. overlap')
                
                #plt.plot(alpha,overlap_changed,'o-',label=str(cloud_numb[m])) 
                #plt.title('overlap vs. normalized height ')
                #plt.xlim([1,0.2])
                
                
                
                
                plt.legend(loc='upper right', bbox_to_anchor=(1.40, 1.0), fontsize='small',ncol=2) 
                
                
                
                
                
                
                
                m = m+1
            
            #mu1=np.mean(Total_dr)
            #sigma1=np.std(Total_dr)
            #Check1=[Total_dr >= mu1 + 1.96 *sigma1]
            
            overlap_c=volume_c/(height_c*projected_area_c)










"""
corr_spear, p_spear = spearmanr(overlap_changed_matrix.flatten(), alpha_matrix.flatten()) #spearman correlation measures monoticity of overlap and normalized ht
spearman=[corr_spear, p_spear]
print('spearman correlation:',corr_spear,', spearman p-value:',p_spear)

avg_dr=np.mean(Total_dr, axis=1) # average dr for each cld
avg_avg_dr=np.mean(Total_dr.flatten()) # avg dr for trackfile, # approx mean of boxplot ?


### calculate interquartile range
#data=Total_dr.flatten()
q25, q75 = percentile(Total_dr.flatten(), 25), percentile(Total_dr.flatten(), 75)
iqr = q75 - q25

### calculate the outlier cutoff
cut_off = iqr * 1.5
lower, upper = q25 - cut_off, q75 + cut_off
#outliers = [x for x in data if x < lower or x > upper]
outliers_great = [x for x in Total_dr.flatten() if x > upper]
### plot box and whisker plot
fig1, ax1 = plt.subplots()
ax1.set_title('Change in Overlap per chop')
ax1.boxplot(Total_dr.flatten())#,notch=True,bootstrap=1000)
### plot histogram
plt.figure()
plt.hist(Total_dr.flatten(), bins='auto')
plt.title('count of change in r vs. change in r per chop')
"""




###################################

end= time.time()
print('Run Time in Seconds:', end-start)

















































