'''
Created on 26 Aug 2015
Geodesy Project to Calculate the offset between predicted ephermeris and post calculated precise ephemeris
@author: craig
'''
from geodesyDec import GeodesyObject


if __name__ == '__main__':


    g=GeodesyObject("../", 85)
    BroadcastDATA = g.readInRinex("docs/brdc2000.15n")
    #_____________GET DATA FOR PRECISE EPHEMERIS_______________________
    #(get precise coordinates AND Range for sat 21 over 24hrs)
    prec_cords,precR =g.readIGS("docs/igs18540.sp3","PG21")
    times=g.getTimes("21")
    for epoch in times:
        #_________GET RANGE AND COORDINATES FOR BROADCAST EPHEMERIS EPOCH
        #(get coordinates for sat 21 at time)
        Br_coords,brdcastR=g.getEarthCenteredCoordinates("21",epoch)

        #_________COMPARE TWO EPHEMERIS FOR EACH EPOCH
        R=g.compare(precR,brdcastR,epoch)
        if epoch == "06 00  0":
            print R



