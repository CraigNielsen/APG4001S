'''
Created on 22 Jul 2015

@author: Craig
'''
from math import factorial, sqrt
import math
import scipy.special

from numpy import cos, sin
from numpy import float_, polynomial
from numpy.polynomial.legendre import legval2d
from scipy.constants.constants import pi

import decimal as DEC
import numpy as np
import sympy as sp


# import fortranformat as ff
DEC.getcontext().prec = 200

class GeodesyObject(object):
    '''
    object for geodesy stuff
    '''


    def __init__(self, dir_,resolution_=160,GM_=3986005.*10**8,a_=6378137.,b_=6356752.3141 ):
        '''
        Constructor
        '''
        self.dir= dir_
        self.resolution=resolution_
        self.GCOEFC1=np.zeros((resolution_,resolution_),dtype=np.dtype(float_))
        self.GCOEFS1=np.zeros((resolution_,resolution_),dtype=np.dtype(float_))
        self.GM=float_(GM_)
        self.a=float_(a_)
        self.b=float_(b_)
        self.RHO=float_(7292115*10**(-11))


        self.gamE=9.7803267715
        self.J2=float_(108236.*10**(-8))


    def  normaliseP(self,n,m,Theta):
        P=self.getP(n,m,Theta)
#         P=self.numpyP(n, m, Theta)
        if m==0:k=1
        else :k=2
        top=k*(2*n+1)*factorial(n-m)
        bottom=factorial(n+m)
         
        Pbar= sqrt((top)/(bottom))*P
        if math.isnan(Pbar):
            print "what"
        return Pbar
    def numpyP(self,n,m,Theta):
        P=scipy.special.lpmn(m,n,Theta)
        print P
        return P 
    
    def getP(self,n,m,Theta):
#         print Theta
        t = sp.symbols('t')
#         print sp.diff((t**2 - 1)**n,t,n)
        funct=1./((2.**n)* factorial(n)) * sp.diff((t**2 - 1)**n,t,n)
        funct2=(1.-t**2)**(m/2.)  * sp.diff(funct,t,m)
        answer=funct2.subs(t,Theta)

#         P=float_(0)
#         t=cos(Theta)
#         fn_= float_(factorial(n))
#         half_= float_((float_(1)  / (2**n *fn_)))
#         half2= float_((1. -t**2)**(m/2))    
#         partOne=    (  half_  * half2 )
#          
#         partTwo=    float_((   (self.getD())**(n+m) *  (t**2-1)**n)/   (self.getD()*t**(n+m)))
#         P=float_(partOne*partTwo)
#         x=legval2d(n,m,Theta)
#         print "still need to check lpmv function here"

#         print answer
        return answer   
            
    def getN(self,r,Latitude,Longitude): 
#         ,r=self.r,Sigma=self.Sig,Lamda=self.Lam,GM=self.GM,a=self.a,
        r=float_(r)
        Lamda=float_(Longitude)
        
        sum=float_(0)
        
        print "CHANGE THIS"
        E = np.sqrt(self.a**2 - self.b**2)
        e = E / self.a
        r = self.a * (np.sqrt(1 - ((((e**2)*(1 - e**2))*((np.sin(Latitude))**2)) /  ((1 - (e**2)*((np.sin(Latitude))**2))))))
        
        
        for n in range(155):
            print n
            if (n<2):continue
            innersum=float_(0)
            for m in range(n):
                C=self.getC(n,m)
                S=self.getS(n,m)
                P=self.normaliseP(n,m,cos(pi/2. - Latitude))
#                 P=self.numpyP(n, m, cos(pi/2. - Latitude))
                l=(cos(float_(m*Lamda)))
#                 print "l is : " + str(l)
                part1=(C* float_(l))
#                 print "part 1: " + str(part1)

                part2=( S * float_(sin(float_(m*Lamda) ) ))
#                 print "part 2 : " + str(part2)
                innersum+=( part1 + part2 )*(P)
#                 print "innersum: " + str(innersum)

            try:
                sum+=float_(((self.a/r)**n)*innersum)
                if math.isnan(sum):
                    print "nooooooo"
            except RuntimeWarning:
                print "nooooooooo"
#             print "sum is : " + str(sum)
        total=  (   (self.GM) /  ( self.getGamma(Latitude)*r ) )* sum
        return total
    
    def getGamma(self,Latitude):
        Lat=pi/2. -Latitude
        g_1 = float_(0.0052790414 *    ((sin(Lat)))**2)
        g_2 = float_(0.0000232718 *    ((sin(Lat)))**4)
        g_3 = float_(0.0000001262 *    ((sin(Lat)))**6)
        g_4 = float_(0.0000000007 *    ((sin(Lat)))**8)
        Gamma = float_(float_(self.gamE )* (1 + g_1 + g_2 + g_3 + g_4)  )
        return Gamma
    
    def getGamA(self):
        
        g= (self.GM / (self.a*self.b))  *   (1- (3./2*self.getM()) - (3./14)*((self.getEdash())**2)*self.getM() )
        
        return g
        
    def getGamB(self):
        g=  (self.GM / (self.a**2.))  *   (1.+ (self.getM()) + (3./7)*((self.getEdash())**2)*self.getM() )
        return g

    def getM(self):
        m= (self.a**2. * self.b * self.RHO**2)/self.GM
        
    def getEdash(self):
        e=(self.a**2.-self.b**2)/self.b**2
        return e
        
    def getC(self,n,m):
        deltaC=float_(self.GCOEFC1[n][m])+self.J2
        return deltaC
    
    def getS(self,n,m):
        S=float_(self.GCOEFS1[n][m])
        return S
    

    
    def getD(self):
        d=1
        return d
    
    def readIn(self,filename):
        np.set_printoptions(precision=22)
        with open(self.dir+filename, "r") as in_file:
            in_line = in_file.readline()

            in_line = in_file.readline()

            while True:
                in_line = in_file.readline()
                if not in_line:
                    break
#                 in_line = in_line[:-1]
                if in_line[5]=="C":
  
                    dat= float_(in_line[24:40])

                    self.GCOEFC1[int(in_line[14:17])][int(in_line[17:22])]=(dat*10**(int(in_line[41:])))
                else:
                    dat= float_(in_line[24:40])
                    self.GCOEFS1[int(in_line[14:17])][int(in_line[17:22])]=(dat*10**int(in_line[41:]))
    
    def writeOut(self,filename):
        fo = open(self.dir+filename, "w")
#         line = ff.FortranRecordWriter('(A1, A1, A20)')
#         line = ff.FortranRecordWriter('(3F15.12)')
        for i in range(161):
            for j in range(161):

                data = self.GCOEFC1[(j)][(i)]
                if (data != 0):
                    fo.write( "GCOEFC1"+ ' {:3.0f} '.format(j) + ' {:3.0f} '.format(i) +"    "+ '%.14g' %data +"\n" )
#                     fo.write(line.write([j,i,data]) + "\n")
#                     fo.write("GCOEFC1"+" "+str(j)+" "+str(i)+" "+("%0.14g" % data)+"\n")
        print "finished writing: "+filename
        fo.close()
        
        
        
        
        
        
        
        
        