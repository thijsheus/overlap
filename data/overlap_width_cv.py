#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 17:03:07 2019

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

#import overlap_calculation

start=time.time()


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
# def 
"""
def overlap(s,h,l):
    #return l / (l +s*h)
    return 1 - (s*h / (l + s*h))
"""
def overlap1(s,h,l,w):
    alpha= 1 / (1+w) + 1
    return 1 - (alpha*0.5*s*h / (l + s*h))

def overlap2(s,h,l):
    #return l / (l +s*h)
    B=1/(1+(h/l))
    return 1 - (B*s*h / (l + B*s*h))

def overlap_inv1(s,h,l):
    #return l / (l +s*h)
    B=h/l
    #return 1 + ( (B*s*h) / l) 
    return  ( (B*s*h) / l)
    
def opt_fun1(X,B):
    #return l / (l +s*h)
    #B=1*h/l
    s,h,l = X
    return  ( (B*s*h) / l) 

def shearpar(X,gs):
    h,l = X
    return gs * (h/l)


####################################
    
# script
filenames=[ricod]
bomexfilenames=[bomextrack360]
ricofilenames=[ricotrack828]
armfilenames=[armtrack522]
conditional_height=0


fractal=  (6*np.sqrt(3)) / (5*np.pi) 

### rico828
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
model_inv=np.load('model_inv_rico2016.npy') 
model_opt=np.load('model_opt_rico2016.npy')
model_opt2=np.load('model_opt2_rico2016.npy')
"""
"""
### bomex 360
height=np.load('height_no_shear_bomex360.npy')
volume=np.load('volume_no_shear_bomex360.npy')
projar=np.load('projarea_no_shear_bomex360.npy')
overlap=np.load('overlap_no_shear_bomex360.npy')
areaz=np.load('area_z_ratio_bomex360.npy')
"""
for file in filenames:
    #zt=file.variables['z'][:]
    zh=file.variables['zh'][:]
    time_t=file.variables['time'][:]
    u=file.variables['u'][:,:]
    v=file.variables['v'][:,:]
    w=file.variables['w'][:,:]

    
            
    if file == ricod:
        for file1 in ricofilenames:
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
            
            nrcloudarray = np.ma.getdata(nrcloud)
            #conditional_height=2000
            index_anvil=np.where(ht > conditional_height)
            index_anvil_size=index_anvil[0].size
            ht_anvil=ht[index_anvil[0]]
            overlap_ratio_anvil=overlap_ratio[index_anvil[0]];
            area_proj_anvil=area_proj[index_anvil[0]]
            print('Clouds that satisfy the given condition: ',index_anvil_size)
            


            index7=index_anvil[0]
            dx=xt[1]-xt[0];dy=yt[1]-yt[0];dz=zt[1]-zt[0];
            gridarea=dx*dy
            gridvol=dx*dy*dz
            nx=xt.size;ny=yt.size;nz=zt.size;
            
            
            cloud_numb=index7 +1 # index is off by 1 as array starts with zero
            
            # selecting cld to plot explicitly
            #cloud_numb=np.array([278,313,351]) #ricotrack2016
            #cloud_numb=np.array([137,371,431]) #ricotrack2016
            #cloud_numb=np.array([70,254,486])   #ricotrack2124
            #cloud_numb=np.array([196,199,242])   #ricotrack2124
            #cloud_numb=np.array([533]) # ricotrack828
            #cloud_numb=np.array([4320, 4526]) # ricotrack2016
            #cloud_numb=np.array([224]) # ricotrack900
            #cloud_numb=np.array([330]) # ricotrack1548
            #cloud_numb=np.array([575]) # ricotrack828
            
            ##################################################### shear
            
            
            
            uz=np.mean(u,axis=0) # avg among diff times
            vz=np.mean(v,axis=0) # avg among diff times
            wz=np.mean(w,axis=0)
            
            UVz= (2*uz*vz)/(uz+vz)
            #UVz= (uz**2 + vz **2)**0.5
            
            wint=interp1d(zh,wz,axis=0)
            wz1=wint(zt)
            
            
            one=np.ones(wz1.size)
            
            duz=np.gradient(uz) # central diff in uz
            dvz=np.gradient(vz)
            
            
            #s_u=duz/np.maximum(wz1,one)
            #s_v=dvz/np.maximum(wz1,one)
            s_u=abs(duz)/np.maximum(wz1,one)
            s_v=abs(dvz)/np.maximum(wz1,one)
            
            
            #Gmean=(s_u*s_v)**(0.5)
            Amean=(s_u+s_v)/2
            #H=Gmean**2/Amean
            Hmean=(2*s_u*s_v) / (s_u + s_v)
            pyth= (s_u**2 + s_v**2)**0.5
            
            shear0=np.zeros(cb.size)
            w_cld=np.zeros(cb.size)
            for i in range(cb.size):
                """
                du=uz[ct[i]]-uz[cb[i]]
                dv=vz[ct[i]]-vz[cb[i]]
                #w1=np.mean(wz1[cb[i]:ct[i]+1])
                su=du #/max(w1,1)
                sv=dv #/max(w1,1)
                #Hmean=(2*su*sv) / max( (su + sv) , 1)
                Amean = (su+sv)/2
                shear0[i]= Amean  #Hmean
                """
                """
                s_u=sum(duz[cb[i]:ct[i]+1]) 
                s_v=sum(dvz[cb[i]:ct[i]+1]) 
                Hmean=(2*s_u*s_v) / (s_u + s_v)
                shear0[i] = Hmean
                """
                shear0[i] = sum(Hmean[cb[i]:ct[i]+1])    # since max(w,1)=1 s_u,s_v = duz, dvz
                
                #shear0[i] = abs(UVz[ct[i]] - UVz[cb[i]])
                
                
                #shear0[i] = np.mean(Hmean[cb[i]:ct[i]+1]) 
                #shear0[i] = Hmean[ct[i]] - Hmean[cb[i]]
                w_cld[i]=np.mean(wz1[cb[i]:ct[i]+1])
            
            """
            bins=2*dz
            
            plt.figure()
            plt.hist2d(overlap_ratio,abs(shear0),bins=bins,cmin=0.5)
            plt.title('shear vs. overlap')
            #plt.ylim([0,10])
            #plt.ylim([min(shear0),max(shear0)])
            colorbar = plt.colorbar()
            colorbar.set_label('counts in bin')
            
            """
            
            
            
            beta_find=( ((1/overlap)-1)*cv**(0.5) ) / (abs(shear0)*ht_anvil  )
            
                
              
            bins=dz
            """
            plt.figure()
            plt.hist2d(overlap1(abs(shear0),ht_anvil,(cv)**(1/2),w_cld),ht_anvil,bins=bins,cmin=0.5)
            plt.title('Height vs. Calculated Overlap')
            plt.xlim([0,1.0])
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            plt.figure()
            plt.hist2d(overlap(abs(shear0),ht_anvil,(cv)**(1/2)),ht_anvil,bins=bins,cmin=0.5)
            plt.title('Height vs. Calculated Overlap')
            plt.xlim([0,1])
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            plt.figure()
            plt.hist2d(overlap_ratio_anvil,ht_anvil,bins=bins,cmin=0.5)
            plt.title('Height vs. Actual Overlap')
            plt.xlim([0,1])
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            """
            """
            plt.figure()
            plt.hist2d(1/overlap1(abs(shear0),ht_anvil,(cv)**(1/2),w_cld),1/overlap_ratio,bins=bins,cmin=2.0)
            plt.plot(np.arange(0,100,1),'k',linewidth=3)
            plt.title('Model Comparison')
            #plt.xlim([0,1.0])
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            plt.xlabel('Calculated Overlap')
            plt.ylabel('Actual Overlap')
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            """
            
            
            shear=   (1/overlap_ratio) - (1/overlap)  
            shear[shear<0]=0
            beta_find=( ((shear))*cv**(0.5) ) / (abs(shear0)*ht_anvil  )
            l=cv**(1/2)
            Model_inv=overlap_inv1(abs(shear0),ht_anvil,l)
            bins=3*dz
            plt.figure()
            ###Note: cmin=5 here while we usually have cmin=2 , or ... cmin=2 and have bins=2*dz 
            plt.hist2d(Model_inv,1/overlap_ratio,bins=[4*dz,1*dz],cmin=3.0) 
            plt.plot(np.arange(0,100,1),'k',linewidth=3)
            plt.title('Model Comparison')
            #plt.xlim([0,1.0])
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            plt.xlim([0,10])
            plt.xlabel('Calculated Overlap')
            plt.ylabel('Actual Overlap')
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            
            
            plt.figure()
            ###Note: cmin=5 here while we usually have cmin=2 , or ... cmin=2 and have bins=2*dz 
            plt.hist2d(shear0,ht,bins=[1*dz,1*dz],cmin=3.0) 
            #plt.plot(np.arange(0,100,1),'k',linewidth=3)
            plt.title('Ht vs. Shear')
            #plt.xlim([0,1.0])
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            #plt.xlim([0,15])
            plt.xlabel('Shear')
            plt.ylabel('Ht')
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            
        
            
        
            
            cmin=3.0
            
            X1=(abs(shear0),ht_anvil,cv**0.5)
            opt1, cov1 = curve_fit(opt_fun1, X1 , shear)
            model_opt1=opt_fun1(X1,*opt1)
            
            plt.figure()
            plt.hist2d(model_opt1, shear, bins=dz , cmin=cmin)
            plt.plot(np.arange(0,100,1),'k',linewidth=3)
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            #plt.xlim([0,20])
            plt.xlabel('Calculated Overlap')
            plt.ylabel('Actual Overlap')
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            
            
            
            """
            
            ####
            #when cmin=3, some (outlier) clds are not shown, this clds have rinv > 5
            ## be careful here , only run once 
            non_out= np.where(model_opt1<5)
            #####
            
            X2=(abs(shear0[non_out]),ht_anvil[non_out],cv[non_out]**0.5)
            opt2, cov2 = curve_fit(opt_fun1, X2 ,1/overlap[non_out])
            model_opt2=opt_fun1((abs(shear0[non_out]),ht_anvil[non_out],cv[non_out]**0.5),*opt2)
            
            #opt2=np.float64([2.4])
            
            plt.figure()
            plt.hist2d(opt_fun1(X2,*opt2), 1/overlap[non_out], bins=dz , cmin=cmin)
            plt.plot(np.arange(0,100,1),'k',linewidth=3)
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            #plt.xlim([0,20])
            plt.xlabel('Calculated Overlap')
            plt.ylabel('Actual Overlap')
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            
            
            shear=opt_fun1(X1,*opt2)
            #r_total =  ( (1/(overlap)) + (1/areaz) + (1/fractal) ) -3  # correction: need to subtract 1 from each 
            #r_total =  ( ((shear)) + (1/areaz) + (1/fractal) ) -3  # correction: need to subtract 1 from each 
            
            
            horizontal= shear
            vertical= 1/overlap
            
            plt.figure()
            plt.hist2d(horizontal, vertical, bins=dz , cmin=cmin)
            plt.plot(np.arange(0,100,1),'k',linewidth=3)
            plt.xlim(left=0)
            plt.ylim(bottom=0)
            #plt.xlim([0,20])
            plt.xlabel('Calculated Overlap')
            plt.ylabel('Actual Overlap')
            colorbar=plt.colorbar(extend='both');
            colorbar.set_label('count in bin')
            plt.clim(0,200)
            plt.xlim([0,5])
            plt.ylim([0,5])
            """
            
end= time.time()
print('Run Time in Seconds:', end-start)

















