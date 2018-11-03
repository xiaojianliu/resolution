# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 12:52:48 2018

@author: xiaojian
"""
from math import radians, cos, sin, atan, sqrt 
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#from scipy import  interpolate
#from datetime import datetime, timedelta
from matplotlib.path import Path
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

def nearest_point(lon, lat, lons, lats, length,x):  #0.3/5==0.06
    '''Find the nearest point to (lon,lat) from (lons,lats),
    return the nearest-point (lon,lat)
    author: Bingwei'''
    #FN='necscoast_worldvec.dat'
    #CL=np.genfromtxt(FN,names=['lon','lat'])
    p = Path.circle((lon,lat),radius=length)
    #plt.figure(figsize=(8,7))
    #numpy.vstack(tup):Stack arrays in sequence vertically
    points = np.vstack((lons.flatten(),lats.flatten())).T  
        
    insidep = []
    #collect the points included in Path.
    for i in np.arange(len(points)):
        if p.contains_point(points[i]):# .contains_point return 0 or 1
            insidep.append(points[i])  
    # if insidep is null, there is no point in the path.
    if len(insidep)<2:
        
        print ('There is no model-point near the given-point.')
        raise Exception()
    #calculate the distance of every points in insidep to (lon,lat)
    distancelist = []
    for i in insidep:
        #print (lon,i[0],'########',lat,i[1])
        ss=haversine(lon, lat, i[0], i[1]) 
        #plt.scatter(i[0],i[1],color='blue',s=1)
        #print (ss)
        distancelist.append(ss)
    mindex = np.argmin(distancelist)
    distancelist[mindex]=100000000
    mindex1 = np.argmin(distancelist)
    #plt.scatter(lon,lat,color='green',s=1)
    #plt.scatter(insidep[mindex1][0],insidep[mindex1][1],color='red',s=1)
    #plt.axis([-70.75,-70,41.72,42.07])
    #plt.plot(CL['lon'],CL['lat'],color='black',linewidth=1)
    #plt.title('cape cod')
    #plt.savefig('cape cod %s'%(x))
    #plt.show()

        
    return distancelist[mindex1]




#url='''http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/archives/necofs_mb?lonc[0:1:165094],latc[0:1:165094]'''
url='''http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?lonc[0:1:90414],latc[0:1:90414]'''

ds = Dataset(url,'r').variables   # netCDF4 version
lonc=ds['lonc'][:]
latc=ds['latc'][:]   
#latitude 41.72,42.07
#longitude -70,-69.8
lon=[]
lat=[]
for a in np.arange(len(lonc)):
    if lonc[a]>-70.75 and lonc[a]<-70 and latc[a]>41.72 and latc[a]<42.07:
        lon.append(lonc[a])
        lat.append(latc[a])

lon=np.array(lon)
lat=np.array(lat)
FN='necscoast_worldvec.dat'
CL=np.genfromtxt(FN,names=['lon','lat'])
plt.figure(figsize=(8,7))
plt.scatter(lon,lat,s=0.5)
plt.axis([-70.75,-70,41.72,42.07])
plt.plot(CL['lon'],CL['lat'],color='black',linewidth=1)
plt.title('cape cod')
plt.savefig('cape cod')
plt.show()

dis=[]
for a in np.arange(len(lon)):
    print (a)
    
    d=nearest_point(lon[a], lat[a], lon, lat, 0.05,a)
    dis.append(d)

meandis=np.mean(dis)
print ('mean distance %s km'%(meandis))
print ('min distance %s km'%(min(dis)))
print ('max distance %s km'%(max(dis)))
