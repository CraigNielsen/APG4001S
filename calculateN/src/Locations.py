'''
Created on 29 Jul 2015

@author: frgcra003
'''
from math import floor

from numpy import radians,pi
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
    
    def addPlace(self,name,lat_,long_,ellipsoidal_height):
        '''NOTE__________INPUT IN [D,M,S] for lat and long'''
        '''______________________________________________'''
        latitude=radians(lat_[0]+lat_[1]/60.+lat_[2]/3600.)
        longitude=radians(long_[0]+long_[1]/60.+long_[2]/3600.)
        r=ellipsoidal_height
        self.place[name]=Place(latitude,longitude,r)
        
    def outputInDegrees(self,name_):
        print name_ + "\n[Lat,Long,r(ellipsoidal)] in [D,M,S]:"
        lat= self.rad2dms(self.place[name_][0])
        long= self.rad2dms(self.place[name_][1])
        r= self.place[name_][2]
        
        return [lat,long,r]
    
    def getRads_LocationObjectFor(self,name_):
#         print name_ + "\n[Lat,Long,r(ellipsoidal)] in Radians:"
        return  self.place[name_]
        
    def rad2dms(self,rad):
        dec=rad*180/pi
        
        deg = int(floor(dec))
        fmin = (dec - deg)* 60.0
        min = int(floor(fmin))
        sec = round(((fmin - min) * 60.0),5)
        return [deg,min,sec]