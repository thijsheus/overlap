#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:47:11 2019

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
bomextrack234 = Dataset('/data/bomex/l.0023400.track.nc','r')
bomextrack252 = Dataset('/data/bomex/l.0025200.track.nc','r')
bomextrack270 = Dataset('/data/bomex/l.0027000.track.nc','r')
bomextrack288 = Dataset('/data/bomex/l.0028800.track.nc','r')
bomextrack306 = Dataset('/data/bomex/l.0030600.track.nc','r')
bomextrack324 = Dataset('/data/bomex/l.0032400.track.nc','r')
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



lassotrack306 = Dataset('/data/lasso/sims/20160611_micro/l.0030600.track.nc','r')


filenames=[bomexd, ricod, armd]

bomexfilenames=[bomextrack18, bomextrack36, bomextrack54, bomextrack72, bomextrack90, bomextrack108, bomextrack126, bomextrack144, bomextrack162, bomextrack180, bomextrack198, bomextrack216]
ricofilenames=[ricotrack36, ricotrack72, ricotrack108, ricotrack144, ricotrack180, ricotrack216, ricotrack252, ricotrack288, ricotrack324, ricotrack360, ricotrack396]
armfilenames=[armtrack108, armtrack126, armtrack144, armtrack162, armtrack180, armtrack198, armtrack216, armtrack234, armtrack252, armtrack270, armtrack288]

###########################################################################

####################################
    
# script
filenames=[bomexd]
ricofilenames=[ricotrack2016]
bomexfilenames=[bomextrack216,bomextrack234,bomextrack252,bomextrack270,bomextrack288,bomextrack306,bomextrack324,bomextrack342,bomextrack360]
armfilenames=[armtrack522]
lassofilenames=[lassotrack306]
conditional_height=0
epsilon1=50    #arbitrarily set, most outliers are > 200, most nonoutliers clds are about 0

conditional_height=0

file1_numb=-1

begin=5
### accessing multiple datasets easily
Afilenames = sorted(glob.glob('/data/bomex/*.track.nc'))
Bfilenames = Afilenames[begin:]
#Afilenames = sorted(glob.glob('/data/rico/*.track.nc'))
#Bfilenames = Afilenames[7:]
#Bfilenames = Afilenames[22:23]

file1_numb=begin-1

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



    index1=np.where(ht > conditional_height) # indices of where condition holds true in ht vector
    index1_size=index1[0].size
    ht=ht[index1[0]]   # taking the ht values according to indices above
    overlap_ratio=overlap_ratio[index1[0]];
    area_proj=area_proj[index1[0]]
    print('Clouds that satisfy the given condition: ',index1_size)
    print('conditional height is',conditional_height,'meters')


    cloud_numb=index1[0] +1 # index is off by 1 as array starts with zero
            
    height_c=np.zeros(cloud_numb.size)
    overlap_c=np.zeros(cloud_numb.size)
    projected_area_c=np.zeros(cloud_numb.size)
    volume_c=np.zeros(cloud_numb.size)
    m=0;
    
    alpha=np.arange(1,0,-1) # just 1
    #alpha=np.arange(1,0.35,-0.1) # fraction
    #alpha=np.array([0.7,0.5,0.3,0.1])
    
    #Dist=np.zeros((cloud_numb.size , alpha.size -1))
    #Dcurve=np.zeros((cloud_numb.size,alpha.size))  # multi dimensional array
    #Dcurve2=np.zeros((cloud_numb.size,alpha.size))  # multi dimensional array
    #Dcurve3=np.zeros((cloud_numb.size,alpha.size))  # multi dimensional array
    
    #Total_dr=np.zeros((cloud_numb.size,alpha.size -1))  # multi dimensional array
    #mu=np.zeros(cloud_numb.size)
    #sigma=np.zeros(cloud_numb.size)
    #Keeptrack_tdr=np.zeros((cloud_numb.size,3))  #make abitrarily big
    #alpha_matrix=np.zeros((cloud_numb.size,alpha.size))
    overlap_changed_matrix=np.zeros((cloud_numb.size,alpha.size))
    area_z_ratio=np.zeros(cloud_numb.size)
    
    for e1 in cloud_numb:  # may take a while to run 
        
        
        location=np.where(nrcloudarray == e1)    # location of all cells with cloud # e1
        layers=(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z 
        #layers=height/dz
        ###
        
        base=np.amin(location[0]);top=np.amax(location[0]);
        #levels=range(base,top+1,1)
        
        x_unique=np.unique(location[2])
        y_unique=np.unique(location[1])
        xtest1=abs(np.mean(x_unique)-np.median(x_unique))  # median is resistant while mean is not
        ytest1=abs(np.mean(y_unique)-np.median(y_unique))
 
        #xtest1=abs(np.mean(location[2])-np.median(location[2]))
        #ytest1=abs(np.mean(location[1])-np.median(location[1]))
        #x_small=np.extract(location[2] < 512,location[2])
        #x_large=np.extract(x_unique > 512,x_unique)
        #y_small=np.extract(y_unique < 512,y_unique)
        #y_large=np.extract(y_unique > 512,y_unique)
        
        if xtest1 > epsilon1:
            splitx=np.where(location[2] < 512)
            location[2][splitx] = location[2][splitx] + nx
            
        if ytest1 > epsilon1:
            splity=np.where(location[1] < 512)
            location[1][splity] = location[1][splity] + ny
        
        
        #sys.exit("Error message")
        

        location_matrix=np.zeros((location[0].size,3))
        location_matrix[:,0]=location[0] # z coor
        location_matrix[:,1]=location[1] # y coor
        location_matrix[:,2]=location[2] # x coor
        COM=np.zeros((top-base+1,3))
        newlocation_matrix=np.zeros(location_matrix.shape)
        cross_area_z=np.zeros(int(layers))
        #u, indices = np.unique(location[0], return_index=True)
        #COM[:,0]=u
        
        
        
        
        
        for z in np.arange(base,top+1,1):
            findz=np.where(location[0] == z)
            cross_area_z[z-base]=dx*dy*len(findz[0])
            # height / dz = layers
            
            
            #COM[z-base,:]=[z,sum(location[1][findz[0]])/findz[0].size, sum(location[2][findz[0]])/findz[0].size]
            COM[z-base,:]=[z,np.mean(location[1][findz[0]]), np.mean(location[2][findz[0]])] #find center of mass w/ arithmetic mean
            ### not quite (close)
            """
            if z<top:
                COM[z-base,1] = sum(location[1][indices[z-base]: indices[z-base+1]-1]) / len(location[1][indices[z-base]: indices[z-base+1]-1])
                COM[z-base,2] = sum(location[2][indices[z-base]: indices[z-base+1]-1]) / len(location[2][indices[z-base]: indices[z-base+1]-1])
            else:
                tempz=np.where(location[0] == z)
                COM[z-base,:]=[z,sum(location[1][tempz[0]])/len(tempz[0]), sum(location[2][tempz[0]])/len(tempz[0])]
        
            """
        area_z_ratio[m]= np.mean(cross_area_z)/np.amax(cross_area_z)    
        center_choice=[np.mean(COM[:,1]), np.mean(COM[:,2])] 
        movement=np.subtract(center_choice , COM[:,1:])
        for z in np.arange(base,top+1,1):
            findz=np.where(location[0] == z)
            newlocation_matrix[findz,1:]=location_matrix[findz,1:] + movement[z-base,:]
        
        newlocation_matrix[:,0]=location_matrix[:,0]
        
        location[0][:]=newlocation_matrix[:,0]   # chang e to calculate overlap w/ code below
        location[1][:]=newlocation_matrix[:,1]
        location[2][:]=newlocation_matrix[:,2]
        
        ########
        
        height_c[m]=dz*(np.amax(location[0]) - np.amin(location[0]) + 1)   # height along z axis
        #projected_area_c[m]=dx*dy*len(np.unique(location[2]+location[1]*nx)) #len of unique # of cells that get proj
        partial_heights=alpha*height_c[m] # fraction of cld height
        #height_diff=height_c[m]-partial_heights
        
        overlap_changed=np.zeros(partial_heights.size)
        projected_area_changed=np.zeros(partial_heights.size)
        vol_changed=np.zeros(partial_heights.size)
        part_loc= -1
        for partial in partial_heights:
            
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
            for j in range(celltop,location[0].size):
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
                """
                 ### plotting in cld in 3d
        
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.scatter3D(dx*location[2][0:celltop+1],dy*location[1][0:celltop+1], dz*location[0][0:celltop+1], c=dz*location[0][0:celltop+1] ,cmap='Paired'); 
                #cmap= single color:copper, cool, winter, multicolor:  Dark2, Paired
                #ax.plot_trisurf(location[2], location[1], location[0], cmap='viridis', edgecolor='none');
                plt.title(str(cloud_numb[m]))
                plt.show()
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
        """
        ## calculates change in r of a cld
        dr=np.zeros(alpha.size -1) 
        for i in range(alpha.size -1 ) :
            dr[i]=  overlap_changed[i+1] - overlap_changed[i]    # having abs lowers the high pt of the interval
            Total_dr[m,i]=dr[i]
            
        """
        """
        alpha_matrix[m,:]=alpha
        overlap_changed_matrix[m,:]=overlap_changed
        Keeptrack_tdr[kt,:]=[file_numb,file1_numb,e1]
        kt=kt+1
        """
        """
        mu[m]=np.mean(dr)
        sigma[m]=np.std(dr)
        
        ### prints index of matrix where change in r is outside the confidence interval
        for i in range(alpha.size -1 ):
            if abs(Total_dr[m,i]) >= mu[m]+1.96*sigma[m]:
                print('index',(m,i))  #remeber indices start at zero
        """
        ###putting alpha = partial_heights / height_c[m] b/c partial_heights = alpha * height_c[m]
        #plt.plot(overlap_changed,alpha,'o-',label=str(cloud_numb[m])) 
        #plt.title('normalized height vs. overlap')
        
        #plt.plot(alpha,overlap_changed,'o-',label=str(cloud_numb[m])) 
        #plt.title('overlap vs. normalized height ')
        #plt.xlim([1,0.2])
        
        
        
        #plt.legend(loc='upper right', bbox_to_anchor=(1.40, 1.0), fontsize='small',ncol=2) 
        
        
        
        
        
        
        
        m = m+1
    
    #mu1=np.mean(Total_dr)
    #sigma1=np.std(Total_dr)
    #Check1=[Total_dr >= mu1 + 1.96 *sigma1]
    
    overlap_c=volume_c/(height_c*projected_area_c)
    
    ################################
    ###
    """


    step=1800 # 30*60
    time_array= np.arange(10800,36000+step,step)
    npyfilespath ="/home/anthonys/Documents/"
    case='bomex'
    
    
    name1=overlap_no_shear_
    name2=volume_no_shear_
    name3=projarea_no_shear_
    name4=height_no_shear_
    name5=aera_z_ratio_
    
    #step=3600 # 60*60
    #time_array= np.arange(10800,216000+step,step)
    #npyfilespath ="/home/anthonys/Documents/"
    #case=rico
    
    ### 
    #np.save(npyfilespath+ name1  + case +str(time_array[file1_numb])+'.npy',overlap_c)
    #np.save(npyfilespath+ name2  + case +str(time_array[file1_numb])+'.npy',volume_c)
    #np.save(npyfilespath+ name3  + case +str(time_array[file1_numb])+'.npy',projected_area_c)
    #np.save(npyfilespath+ name4  + case +str(time_array[file1_numb])+'.npy',height_c)
    #np.save(npyfilespath+ name5  + case +str(time_array[file1_numb])+'.npy',area_z_ratio_c)
    
    
    #np.save(npyfilespath+str(time_array[20])+'.npy',area_z_ratio)
    ### need to fix  saving
    """

    ###########################


plt.plot(overlap_c,overlap_ratio,'o')
plt.plot(overlap_c,overlap_c)
#plt.title('Actual Overlap vs. Adjusted Overlap')
plt.title('Actual Overlap vs. Overlap with no shear')












###################################

end= time.time()
print('Run Time in Seconds:', end-start)





























