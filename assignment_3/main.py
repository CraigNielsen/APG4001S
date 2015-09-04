'''
Created on 26 Aug 2015

@author: craig
'''
from geodesyDec import GeodesyObject


if __name__ == '__main__':
    g=GeodesyObject("../", 85)
    dict = g.readInRinex("docs/brdc2000.15n")
    for name in dict.keys():
        if name =="21":
            for time in dict[name].keys():
                print time
        
    # precD =g.readIGS("docs/igs18540.sp3")
    # for sat in precD.keys():
    #     print precD[sat][" 0  0"]["x"]
    #     print precD[sat][" 0  0"]["y"]
    #     print precD[sat][" 0  0"]["z"]