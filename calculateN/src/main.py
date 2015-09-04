'''
Created on 22 Jul 2015

@author: Craig
'''
from numpy import float_
from GeodesyDecimal import GeodesyObject
from numpy import float_,polynomial,cos,pi,sin
from numpy.polynomial.legendre import *
from Locations import Locations_
from scipy import * 
import scipy.special
from place import Place
from scipy.special import eval_legendre
import numpy as np

import numpy

def coefficients(order):
    for i in range(1, order):
        base = np.zeros(order)
        base[i] = 1
        yield base
    
def a(n, m):
    return 1.0*(2*n+1) / ((n*(n+1))**m)
    
def g(const_dist, m, order):
    legendres = [np.polynomial.legendre.Legendre(n) for n in coefficients(order)]
    terms = [a(n+1,m)*legendres[n](const_dist) for n,_ in enumerate(legendres)]
    return sum(terms)
 
if __name__ == '__main__':

    n=0
    m=0
    Theta=cos(5)

    loc=Locations_()
    loc.addPlace("Hermanus", [34,25,28.6671], [19,13,23.0264], 63.048)
#     loc.addPlace("Pretoria",[25,43,55.2935],[28,16,57.4873],1387.339)
#     loc.addPlace("Richards Bay",[28,47,43.9616],[32,4,42.1896],31.752)
#     loc.addPlace("Thohoyandou",[23,4,47.6714],[30,23,2.4297],630.217)
#     loc.addPlace("Ulundi",[28,17,35.2196],[31,25,15.3309],607.947)
# 
    gO=GeodesyObject("../", 161)
#     gO.readIn("fileIO.txt")
# # #     gO.writeOut("food.txt")
#     hermanus=loc.getRads_LocationObjectFor("Hermanus")
#     print (gO.getN(hermanus.r,  hermanus.latitude, hermanus.longitude))
    
    print eval_legendre(5, cos(pi/8.))
    print gO.normaliseP(5, 5, cos(pi/8.))
    print g((pi),1,1)
#     print gO.numpyP(5,2,cos(pi/8.))


    
#     value=pi/8
#     m=1
#     n=1
#     
#     print gO.normaliseP(n, m, cos(value))
# 
#     print -sin(value)
#     a = 6378137.0             
#     b = 6356752.3141  
#     u = sin(((pi)/2)-(arctan(((b/a)**2)*tan(value))))
#     print "u is " + str( u)
#     loc["home"] = [10, 10, 100]
    
#     print loc.outputInDegrees("Hermanus")
#     print loc.getRadians("Hermanus")