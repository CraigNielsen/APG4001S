'''
Created on 3 Aug 2015

@author: Craig
'''
from numpy import sqrt, sin


class Place(object):
    '''
    classdocs
    '''


    def __init__(self, latitude_,longitude_,ellipsoidal_height_,list=[0,0,0],name_="temp"):
        '''
        a capsule for a place (latitude, longitude and Radius from earth)
        '''
        self.latitude=latitude_
        self.longitude=longitude_
        self.r=ellipsoidal_height_
        self.x=list[0]
        self.y=list[1]
        self.z=list[2]
        self.name=name_
        
    def getR(self):
        a=6378137.0
        b=6356752.3141
        E=sqrt(a**2-b**2)
        e2=(E/a)**2
        lat=self.latitude
        
#         R = a * (sqrt(1 -   (   (  ( (e2)*(1 - e2)) *    ((sin(lat))**2)) / ((1 - (e2)*((sin(lat))**2))))    ))
        return sqrt(self.x**2 + self.y**2 + self.z**2)
        