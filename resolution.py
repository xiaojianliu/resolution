# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 12:52:48 2018

@author: xiaojian
"""
from math import radians, cos, sin, atan, sqrt 
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import  interpolate
from datetime import datetime, timedelta
def haversine(lon1, lat1, lon2, lat2): 
    """ 
    Calculate the great circle distance between two points  
    on the earth (specified in decimal degrees) 
    """   
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])  
    #print 34
    dlon = lon2 - lon1   
    dlat = lat2 - lat1   
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2  
    c = 2 * atan(sqrt(a)/sqrt(1-a))   
    r = 6371 
    d=c * r
    #print type(d)
    return d

url='''http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?lonc[0:1:90414],latc[0:1:90414]'''

ds = Dataset(url,'r').variables   # netCDF4 version
lonc=ds['lonc'][:]
latc=ds['latc'][:]   
#latitude 41.5,42.23
#longitude -70.75,-69.8
lon=[]
lat=[]
for a in np.arange(len(lonc)):
    if lonc[a]>-70.75 and lonc[a]<-69.8 and latc[a]>41.5 and latc[a]<42.23:
        lon.append(lonc[a])
        lat.append(latc[a])
plt.figure(figsize=(7,7))
plt.scatter(lon,lat,s=1)
lon1=-70.4
lat1=42
d=[]
for a in np.arange(len(lon)):
    d.append((lon[a]-lon1)*(lon[a]-lon1)+(lat[a]-lat1)*(lat[a]-lat1))
index1=np.argmin(d)
d[index1]=100000000
index2=np.argmin(d)
plt.scatter(lon[index1],lat[index1],color='green')
plt.scatter(lon[index2],lat[index2],color='red')

lon2=lon[index1]
lat2=lat[index1]

lon3=lon[index2]
lat3=lon[index2]
dd=haversine(lon2,lat2,lon3,lat3)
print ('resolution',dd)

