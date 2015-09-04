'''
Created on 29 Jul 2015

@author: frgcra003
'''
from math import floor
from numpy import radians, pi, sqrt

from place import Place


class Locations_(dict):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.place={}

        self.x=0.
        self.y=0.
        self.z=0.
    
    def addPlace(self,name,lat_,long_,ellipsoidal_height,list_=[0,0,0]):
        '''NOTE__________INPUT IN [D,M,S] for lat and long'''
        '''______________________________________________'''
        latitude=radians(lat_[0]+lat_[1]/60.+lat_[2]/3600.)
        longitude=radians(long_[0]+long_[1]/60.+long_[2]/3600.)
        r=ellipsoidal_height
        list=list_
        self.place[name]=Place(latitude,longitude,r,list,name)
        
    def outputInDegrees(self,name_):
        print name_ + "\n[Lat,Long,r(ellipsoidal)] in [D,M,S]:"
        lat= self.rad2dms(self.place[name_][0])
        long= self.rad2dms(self.place[name_][1])
        r= self.place[name_][2]
        
        return [lat,long,r]
    
    def getRads_LocationObjectFor(self,name_):
        self.name=name_
#         print name_ + "\n[Lat,Long,r(ellipsoidal)] in Radians:"
        return  self.place[name_]
        
    def rad2dms(self,rad):
        dec=rad*180/pi
        
        deg = int(floor(dec))
        fmin = (dec - deg)* 60.0
        min = int(floor(fmin))
        sec = round(((fmin - min) * 60.0),5)
        return [deg,min,sec]
    
    def addCartesian(self,xyzlist):
        self.x=xyzlist[0]
        self.y=xyzlist[1]
        self.z=xyzlist[2]
        
        
