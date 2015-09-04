'''
Created on 22 Jul 2015

@author: Craig
'''
from math import factorial
from numpy import cos, sin
from scipy import float_

import decimal as DEC
import fortranformat as ff
import numpy as np


DEC.getcontext().prec = 100

class GeodesyObject(object):
    '''
    object for geodesy stuff
    '''


    def __init__(self, dir_,resolution_=160,GM_=39860004.418*10**8,a_=6378137.,b_=6356752. ):
        '''
        Constructor
        '''
        self.dir= dir_
        self.resolution=DEC.Decimal(resolution_)
        self.GCOEFC1=np.zeros((resolution_,resolution_),dtype=np.dtype(DEC.Decimal))
        self.GCOEFS1=np.zeros((resolution_,resolution_),dtype=np.dtype(DEC.Decimal))
        self.GM=DEC.Decimal(GM_)
        self.a=DEC.Decimal(a_)
        self.b=DEC.Decimal(b_)
        self.RHO=DEC.Decimal(7292115*10**(-11))
        self.r=DEC.Decimal(0)
        self.Lam=DEC.Decimal(0)


        
    def getG(self,r,Theta,Lamda): 
#         ,r=self.r,Sigma=self.Sig,Lamda=self.Lam,GM=self.GM,a=self.a,
        self.r=r
        self.Lam=Lamda
        
        sum=DEC.Decimal(0)
        innersum=DEC.Decimal(0)
        for n in range(self.resolution):
            if (n<2):continue
            for m in range(n):
                C=self.getC(n,m)
                S=self.getS(n,m)
                P=self.getP(n,m,Theta)
                part1=(C* DEC.Decimal(cos(m*Lamda)))
                part2=( S * DEC.Decimal(sin(m*Lamda)  ))
                innersum=( part1 + part2 )*(P)
                
            sum+=((self.a/r)**n)*innersum
        total=  (   (self.GM) /  (self.getGamma*r) )* sum
        return total
    
    def getGamma(self):
        Gamma=0
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
        C=DEC.Decimal(self.GCOEFC1[n][m])
        return C
    
    def getS(self,n,m):
        S=DEC.Decimal(self.GCOEFS1[n][m])
        return S
    
    def getP(self,n,m,Theta):
        P=DEC.Decimal(0)
        t=cos(Theta)
        fn_= DEC.Decimal(factorial(n))
        half_= DEC.Decimal((DEC.Decimal(1)  / (2**n *fn_)))
        half2= DEC.Decimal((1. -t**2)**(m/2))    
        partOne=    (  half_  * half2 )
        
        partTwo=    DEC.Decimal((   (self.getD())**(n+m) *  (t**2-1)**n)/   (self.getD()*t**(n+m)))
        P=DEC.Decimal(partOne*partTwo)
        return P
    
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
        line = ff.FortranRecordWriter('(3F15.12)')
        for i in range(161):
            for j in range(161):

                data = self.GCOEFC1[(j)][(i)]
                if (data != 0):
                    fo.write( "GCOEFC1"+ ' {:3.0f} '.format(j) + ' {:3.0f} '.format(i) +"    "+ '%.14g' %data +"\n" )
#                     fo.write(line.write([j,i,data]) + "\n")
#                     fo.write("GCOEFC1"+" "+str(j)+" "+str(i)+" "+("%0.14g" % data)+"\n")
        print "finished writing: "+filename
        fo.close()
        
        
        
        
        
        
        
        
        