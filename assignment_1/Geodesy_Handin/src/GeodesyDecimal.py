'''
Created on 22 Jul 2015

@author: Craig
'''
from math import factorial
import math
from numpy import *

import numpy as np
import sympy as sp


# import fortranformat as ff
class GeodesyObject(object):
    '''
    object for geodesy stuff
    '''
    def __init__(self, dir_,resolution_=161,GM_=3986005.0e8,a_=6378137.,b_=6356752.3141 ):
        self.dir= dir_
        self.resolution=resolution_
        self.GCOEFC1=np.zeros((161,161),dtype=np.dtype(float_))
        self.GCOEFS1=np.zeros((161,161),dtype=np.dtype(float_))
        self.GM=float_(GM_)
        self.a=float_(a_)
        self.b=float_(b_)
        self.RHO=float_(7292115*10**(-11))
        self.gamE=9.7803267715
        self.J2 =108263*(10**-8)
        self.J4 =-0.00000237091222
        self.J6 =0.00000000608347
        self.J8 =-0.00000000001427
        self.P={}
        self.P_normalised={}
        
        
    def  normaliseP(self,n,m,Theta):
         
        P=self.getP(n,m,Theta)
 
#         P=self.numpyP(n, m, Theta)
        if m==0:k=1
        else :k=2
         
        top=float_(k*(2*n+1)*math.factorial(n-m))
        bottom=float_(math.factorial(n+m))
          
        Pbar= np.sqrt((top)/(bottom))*P
        if math.isnan(Pbar):
            print "error in NormaliseP"
        return Pbar

    
    def getP(self,n,m,Theta):
        t_ = np.cos(( (np.pi)/2) -       (np.arctan(((self.b/self.a)**2)*np.tan(  Theta))) )
        t = sp.symbols('t')
#         print sp.diff((t**2 - 1)**n,t,n)
        funct=1./((2.**n)* math.factorial(n)) * sp.diff((t**2 - 1)**n,t,n)
        funct2=(1.-t**2)**(float_(m)/2.)  * sp.diff(funct,t,m)
        
        answer=funct2.subs(t,t_)
 
        return answer  
    
    def getPtext(self,name):
        with open(name, "r") as in_file:
            while True:
                in_line = in_file.readline()
                if not in_line:
                    break
#                 in_line = in_line[:-1]

                n=int(in_line[2:5])
                m=int(in_line[9:11])
#                 dat=float_(in_line[12:30])
                dat=float( (math.ceil(float_(in_line[12:30])*1000000000))/1000000000)
                item=int(str(n)+str(m))
#                 print n
#                 print m                
#                 print dat

                self.P_normalised[item]=dat
                
    def getN(self,location_,PrecisionOfN): 

        self.resolution=PrecisionOfN
        Lamda=float_(location_.longitude)           #LAMDA IS THE LONGITUDE IN RADIANS
#         self.getNormalisedP( location_.latitude)   #changed to Colatitude inside here QUESTION
        self.getPtext(location_.name+".txt")
        
        sum=float_(0)
        totalN={}
        r=float_(location_.getR())

        for n in range(2,self.resolution ):
            innersum=float_(0)
            for m in range(0,n+1):
                item = int(str(n) + str(m))
                C=self.getC(n,m)
                C=float_(round(C,12))
                S=self.getS(n,m)
                S=float_(round(S,12))
                if S.round(12)==0 or C.round(12)==0:
                    continue

                P=self.P_normalised[item]
#                 P=self.normaliseP(n, m, location_.latitude)

                innersum+=(((C*cos(m*Lamda))+(S*sin(m*Lamda)))*P)
            sum+=float_((math.pow(self.a/r,n))*innersum)
            
        total=  (   (self.GM) / ( ( self.getGamma(location_.latitude)*r )) )* sum
        return total
#__________________________________________________________________________
#_______________P Normalised Calculation 1 ______________________________________________
    
    def getNormalisedP(self,lat_):

        t = np.cos(( (np.pi)/2) -          (  lat_)) #using colatitude of geocentric latitude
#         u = np.sin((    (np.pi)/2) -    (np.arctan(((self.b/self.a)**2)*np.tan( lat_))) )



        for n in range(0,self.resolution ):    
            N=0    
            for m in range(0,n+1):
                r = (n - m) / 2
                Sum = 0.0   
                for k in range (0,r+1):
                    left = (-1)**k
                    top= math.factorial(2*n - 2*k)
                    bottom= math.factorial(k)*math.factorial(n-k)*math.factorial(n-m-2*k)
                    right= t**(n-m-2*k)
                    Sum += left*    (float_(top)/float_(bottom)) *right

                Pnm_=((2**(-n)*(1-t**2)**(m/2.)  * Sum))
                if (m==0):
                    k=1
                else:             
                    k=2
#__________________________________________________________________________
#_______________Normalization______________________________________________                    
                               
                p1 = float_(k * (2.0*n+1.0))
                p2 = float_(factorial(n-m))
                p3 = float_(factorial(n+m))
                p4 = p2/p3
                PnmNorm = sqrt(p1*p4) * Pnm_
                item = int(str(n) + str(m))
                self.P_normalised[item] = PnmNorm    
      
            

    def getGamma(self,Latitude):
        Lat=pi/2. - Latitude
        k=0.001931851353
        top=1+ k*(sin(Lat)**2)
        bottom=sqrt(1-0.00669438002290*(np.power(sin(Lat),2)))
        gamma=self.gamE* (top/bottom)
        return gamma
    

        
    def getC(self,n,m):
    
        if (n==2 and m==0):
            return float_(self.GCOEFC1[n][m])-self.J2
        elif (n==4 and m==0):
            return float_(self.GCOEFC1[n][m])-self.J4
        elif (n==6 and m==0):
            return float_(self.GCOEFC1[n][m])-self.J6
        elif (n==8 and m==0):
            return float_(self.GCOEFC1[n][m])-self.J8
        else: return float_(self.GCOEFC1[n][m])
    
    def getS(self,n,m):
        S=float_(self.GCOEFS1[n][m])
        return S
    
    def readIn(self,filename):
#         np.set_printoptions(precision=22)
        with open(self.dir+filename, "r") as in_file:
            in_line = in_file.readline()
            in_line = in_file.readline()
            while True:
                in_line = in_file.readline()
                if not in_line:
                    break
#                 in_line = in_line[:-1]
                if in_line[5]=="C":
  
                    dat= str(in_line[24:])
                    dat=(dat.replace("D", "e"))
                    dat=float_(dat)
                    self.GCOEFC1[int(in_line[14:17])][int(in_line[17:22])]=(dat)
                elif in_line[5]=="S":
                    dat= str(in_line[24:])
                    dat=(dat.replace("D", "e"))
                    dat=float_(dat)
                    self.GCOEFS1[int(in_line[14:17])][int(in_line[17:22])]=(dat)
                else:print "error"
    def writeOut(self,filename):
        fo = open(self.dir+filename, "w")
        for i in range(161):
            for j in range(161):
                data = self.GCOEFC1[(j)][(i)]
                sdata=self.getS(j,i)
                if (data != 0):
                    fo.write( "GCOEFC1"+ ' {:3.0f} '.format(j) + ' {:3.0f} '.format(i) +"    "+ '{:22.16g}'.format(sdata) +"   "+ '%.14g' %data +"\n" )
        print "finished writing: "+filename
        fo.close()
        
        
        
    def get_stegunNormaisedP(self,Lat_):
        normLat=( (np.arctan(((self.b/self.a)**2)*np.tan(      Lat_))))
        print "check above pi/2 missing in getHolmes"
        t = np.sin(normLat)  
        u = np.cos(normLat)
        P={}
        P[00] = 1         
        P[10] = t
        P[11] = u
        P[20] = (1.5)*(t**2) - (0.5)
        P[21] = 3* u * t
        P[22] = 3 * (u**2)
        
        
       
        def findP1(n):

            if n==2:
                return P[20]
            elif n==1:
                return P[10]
            elif n==0:
                print "fixed this"
                return P[00]
            else:
                return (2*n+1)*t*findP1(n-1)-n*findP1(n-2)
        
        print findP1(5)

        
        def findPnn(n,m):
            if n==2 and m==0:
                return P[20]
            elif n==2 and m==1:
                return P[21]
            elif n==2 and m==2:
                return P[22]
            elif n==1 and m==0:
                return P[10]
            elif n==1 and m==1:
                return P[11]
            elif n==0 and m==0:
                print "fixed this"
                return P[00]
            
            else:
                return (2*n-1)*u*findPnn(n-1,n-1)
            
        def findPnm(n,m):
            if n==2 and m==0:
                return P[20]
            elif n==2 and m==1:
                return P[21]
            elif n==2 and m==2:
                return P[22]
            elif n==1 and m==0:
                return P[10]
            elif n==1 and m==1:
                return P[11]
            elif n==0 and m==0:
                print "fixed this"
                return P[00]
            else:
                return findPnm(n-2,m)+(2*n-1)*u*findPnm(n-1,m-1)
            
        for n in range (2,self.resolution):
            for m in range (0,n+1):
                 
                item=int(str(n)+str(m))
                if m==0:k=1
                elif m!=0:k=2
                norm=math.sqrt((k*(2*n +1)*(math.factorial(n-m)/math.factorial(n+m))))
                if norm == 0 : print "zero"
                 
                if m==0:
                    item=int(str(n+1)+str(m))
                    P[item]=norm * findP1(n)
                elif n==m:
                    P[item]=norm * findPnn(n, m)
                else:
                    P[item]=norm * findPnm(n, m)
        self.P=P
#                 i1  =   int(str(n)+str(0)) 
#                 i2  =   int(str(n-1)+str(0))   
#                 P[item]=(2*n+1)* t * P[i1] - n*P[i2]
#                 
#                 item=int(str(n)+str(n))
#                 i1  =   int(str(n-1)+str(n-1))   
#                 P[item]=(2*n-1)* u * P[i1]
#                 
#                 if not (m==0) and not (n-2<m):
#                     item=int(str(n)+str(m))
#                     i1  =   int(str(n-2)+str(m))
#                     i2  =   int(str(n-1)+str(m-1))   
#                     P[item]=P[i1] + (2*n-1)* u * P[i2]
#         self.P=P
        
        
    def get_HolmesNormalisedP(self,lat_):
        normLat=(  (np.arctan(((self.b/self.a)**2)*np.tan( lat_))))
        print "check above in getHolmes"
        t = np.sin(normLat)  
        u = np.cos(normLat)
        P={}
        P[00] = 1         
        P[10] = t
        P[11] = u
        P[20] = (1.5)*(t**2) - (0.5)
        P[21] = 3* u * t
        P[22] = 3 * (u**2)

        for n in range(0,2):        
            for m in range(0,n+1):
                if m==0:k=1
                if m!=0:k=2
                norm=sqrt((k*(2*n +1)*math.factorial(n-m)/math.factorial(n+m)))
                item=int(str(n)+str(m))
                P[item]=norm*P[item]
                print "p is: " + str( P[item])
        
        
        P[00] = 1         
        P[11] = sqrt(3)*u

        
        def PnonSect(n,m):
            if n==0 and m==0:
                return P[00]
            elif n==1 and m==1:
                return P[11]
            elif n<m:
                return 0
            else:
                anm=sqrt(   float_((2*n-1)*(2*n + 1)) / float_((n-m)*(n+m))  )
                bnm=sqrt(   float_((2*n+1)*(n+m-1)*(n-m-1))  /   float_((n-m)*(n+m)*(2*n-3))    )
                return anm*t*PnonSect(n-1, m)-bnm*PnonSect(n-2, m)
            
        def PSect(n,m):
            if n==0 and m==0:
                return P[00]
            elif n==1 and m==1:
                return P[11]
            elif n<m:
                return 0
            else:
                return u*math.sqrt((2*m + 1)/2*m)*PSect(n-1, m-1)
                
               
        for n in range (2,self.resolution):
            for m in range (0,n+1):
                 
                item=int(str(n)+str(m))
                if n>m:
                    P[item]=PnonSect(n, m)
                elif n==m:
                    P[item]=PSect(n, m)
        
        self.P=P     
        
        
        
        
        