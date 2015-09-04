'''
Created on 22 Jul 2015

@author: Craig
'''

from GeodesyDecimal import GeodesyObject
from Locations import Locations_


if __name__ == '__main__':



    loc=Locations_()
    x=4973168.840
    y=1734085.512
    z=-3585434.051
 
    loc.addPlace("Hermanus", [-34,-25,-28.6671], [19,13,23.0264], 63.776,list_=[x,y,z])
    loc.addPlace("Pretoria",[-25,-43,-55.2935],[28,16,57.4873],1387.339,[5064032.237,2724721.031,-2752950.762])
    loc.addPlace("Richards Bay",[-28,-47,-43.9616],[32,4,42.1896],31.752,[4739765.776,2970758.460,-3054077.535])
    loc.addPlace("Thohoyandou",[-23,-4,-47.6714],[30,23,2.4297],630.217,[5064840.815,2969624.535,-2485109.939])
    loc.addPlace("Ulundi",[-28,-17,-35.2196],[31,25,15.3309],607.947,[4796680.897,2930311.589,-3005435.714])
    
    gO=GeodesyObject("../")
    gO.readIn("fileIO.txt")

    hermanus=loc.getRads_LocationObjectFor("Hermanus")
    Pretoria=loc.getRads_LocationObjectFor("Pretoria")
    Richards=loc.getRads_LocationObjectFor("Richards Bay")
    Thohoyandou=loc.getRads_LocationObjectFor("Thohoyandou")
    Ulundi=loc.getRads_LocationObjectFor("Ulundi")
#     fo = open("Pretoria.txt", "w")    
    print "hermanus: ", (gO.getN(hermanus,86))
    print "Pretoria: ", (gO.getN(Pretoria,86))
    print "Richards: ", (gO.getN(Richards,86))
    print "Thohoyandou: ", (gO.getN(Thohoyandou,86))
    print "Ulundi: ", (gO.getN(Ulundi,86))

#     for n in range (0,86):
#         print n
#         for m in range(0,n+1):
#             P=gO.normaliseP(n, m, hermanus.latitude)
#             fo.write( "n:{:3} m:{:3} {:3.10g} ".format(n,m,P)+"\n" )
#     fo.close()

#     fo = open("Richards.txt", "w")
#     for n in range (0,86):
#         print n
#         for m in range(0,n+1):
#             P=gO.normaliseP(n, m, hermanus.latitude)
#             fo.write( "n:{:3} m:{:3} {:3.10g} ".format(n,m,P)+"\n" )
#     fo.close()

#     fo = open("Thohoyandou.txt", "w")
#     for n in range (0,86):
#         print n
#         for m in range(0,n+1):
#             P=gO.normaliseP(n, m, hermanus.latitude)
#             fo.write( "n:{:3} m:{:3} {:3.10g} ".format(n,m,P)+"\n" )
#     fo.close()

#     fo = open("Ulundi.txt", "w")
#     for n in range (0,86):
#         print n
#         for m in range(0,n+1):
#             P=gO.normaliseP(n, m, hermanus.latitude)
#             fo.write( "n:{:3} m:{:3} {:3.10g} ".format(n,m,P)+"\n" )
#     fo.close()
#     print (gO.getN(Pretoria,86))
#     print (gO.getN(Richards,86))
#     print (gO.getN(Thohoyandou,86))
#     print (gO.getN(Ulundi,86))


