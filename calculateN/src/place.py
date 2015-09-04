'''
Created on 3 Aug 2015

@author: Craig
'''

class Place(object):
    '''
    classdocs
    '''


    def __init__(self, latitude_,longitude_,ellipsoidal_height_):
        '''
        a capsule for a place (latitude, longitude and Radius from earth)
        '''
        self.latitude=latitude_
        self.longitude=longitude_
        self.r=ellipsoidal_height_
        
        