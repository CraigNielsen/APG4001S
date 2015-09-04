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
        self.J2 =108263.0e-8
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
        
        top=float_(k*(2*n+1)*factorial(n-m))
        bottom=factorial(n+m)
         
        Pbar= sqrt((top)/(bottom))*P
        if math.isnan(Pbar):
            print "error in NormaliseP"
        return Pbar

    
    def getP(self,n,m,Theta):

        t = sp.symbols('t')
#         print sp.diff((t**2 - 1)**n,t,n)
        funct=1./((2.**n)* factorial(n)) * sp.diff((t**2 - 1)**n,t,n)
        funct2=(1.-t**2)**(m/2.)  * sp.diff(funct,t,m)
        answer=funct2.subs(t,Theta)

        return answer  
    def numpyP(self,n,m,Theta):
        P=scipy.special.lpmn(m,n,Theta)

        return P 
    
    def getAllFortranP(self,Theta):
        for n in range(2,self.resolution):        
            for m in range(0,n+1):

                r = (n - m) / 2
                modr = (n-m)%2
                if modr != 0:
                    r = (n - m - 1)/2
                Sum = 0.0   
                
                for k in range (0,r+1):
                    a = (-1)**k
                    b = (n - m - 2*k)
                    c = math.factorial(2*n - 2*k)
                    d = math.factorial(k)
                    f = math.factorial(n - k)
                    g = math.factorial(n - m - 2*k) 
                    h = Theta ** (n - m -2*k)
                    
                    this_loop = a* (c/(d*f*g))*h
        
                    Sum += this_loop
        
                Pnm = (2 ** (-n) * (1 - Theta ** 2) ** (m/2.0)) * Sum
                item = int(str(n) + str(m))
                self.P[item] = Pnm
                
    def getAllFortanPNormalised(self):

        for n in range(2,self.resolution):
            for m in range (0, n+1):
                item = int(str(n) + str(m))
                P = self.P[item]
                if m == 0: k = 1
                else:   k = 2
                p1 = k * (2.0*n+1.0)
                p2 = float_(factorial(n-m))
                p3 = float_(factorial(n+m))
                p4 = p1*p2/p3
        
                self.P_normalised[item] = (sqrt(p4)) * P
        
            
    def getN(self,location_): 


        Lamda=float_(location_.longitude)
        self.getAllFortranP(cos(pi/2. - location_.latitude))
        self.getAllFortanPNormalised()
        sum=float_(0)
        
        print "CHANGE THIS"
        E = np.sqrt(self.a**2 - self.b**2)
        e = E / self.a
        r = self.a * (np.sqrt(1 - ((((e**2)*(1 - e**2))*((np.sin(location_.latitude))**2)) /  ((1 - (e**2)*((np.sin(location_.latitude))**2))))))
        
        r=float_(location_.getR())
        print "r is: " + str(r)
        for n in range(self.resolution):
            if (n<2):continue
            innersum=float_(0)
            for m in range(n+1):
                item = int(str(n) + str(m))
                C=self.getC(n,m)
                S=self.getS(n,m)
                P=self.P_normalised[item]
#                 P=self.normaliseP(n,m,cos(pi/2. - Latitude))

                l=(cos(float_(m*Lamda)))


                part1=(C* float_(l))


                part2=( S * float_(sin(float_(m*Lamda) ) ))

                innersum+=( part1 + part2 )*(P)


            try:
                sum+=float_(((self.a/r)**n)*innersum)
                if math.isnan(sum):
                    print "nooooooo"
            except RuntimeWarning:
                print "nooooooooo"
#             print "sum is : " + str(sum)
        total=  (   (self.GM) /  ( self.getGamma(location_.latitude)*r ) )* sum
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
        if (m==0):
            if (n==2):
                return float_(self.GCOEFC1[n][m])+self.J2 
            elif (n==4):
                return float_(self.GCOEFC1[n][m])+self.J4
            elif (n==6):
                return float_(self.GCOEFC1[n][m])+self.J6
            elif (n==8):
                return float_(self.GCOEFC1[n][m])+self.J8
            
            else:return float_(self.GCOEFC1[n][m])

        else:
            return float_(self.GCOEFC1[n][m])
    
    def getS(self,n,m):
        S=float_(self.GCOEFS1[n][m])
        return S
    

    
    def getD(self):
        d=1
        return d
    def readIGS(self,filename):
        
        with open (filename,"r") as file:
            lines = file.readlines()
            i=0
            dict_={}


            while i < len(lines):
                if lines[i][0]=="*":
                    time=lines[i][14:19]

                    i+=1
                if not lines[i][0]=="P": i+=1 ; continue
                #________________________________________
                
                name=lines[i][:4]
                inDic={}
                timeDic={}
                inDic["x"]=float_(lines[i][4:19])
#                 print lines[i][4:19]
                inDic["y"]=float_(lines[i][20:33])
                inDic["z"]=float_(lines[i][34:48])
                
                
                if name in dict_:
                    dict_[name][time]=inDic 
                else:   
                    timeDic[time]=inDic
                    dict_[name]=timeDic
                i+=1
        return dict_
    
    def readInRinex(self,filename):
        endOfHeader=False
        space=19
        dict_={}
        with open(filename,"r") as file:
            lines=file.readlines()
            i=0
            while i < len(lines):
    
                if lines[i].strip()=="END OF HEADER":
                    endOfHeader=True
                elif not endOfHeader: 
                    i+=1
                if endOfHeader:
                    #for every new point : new name and time
                    if not lines[i][0] == " ":
                        name = lines[i][:2]
                        time = lines[i][12:20]

                        inDic={}
                        timeDic={}
                        inDic["svBias"]=lines[i][22:41]
                        inDic["svDrift"]=lines[i][42:61]
                        inDic["svDRate"]=lines[i][61:81]
                        i+=1
                        
                        inDic["IODE"]=lines[i][3:22]
                        inDic["crs"]=lines[i][22:41]
                        inDic["DelN"]=lines[i][41:60]
                        inDic["MO"]=lines[i][60:81]
                        i+=1
                        inDic["Cuc"]=lines[i][3:22]
                        inDic["e"]=lines[i][22:41]
                        inDic["Cus"]=lines[i][41:60]
                        inDic["sqrtA"]=lines[i][60:81]
                        i+=1
                            
                        inDic["Toe"]=lines[i][3:22]
                        inDic["Cic"]=lines[i][22:41]
                        inDic["OMEGA"]=lines[i][41:60]
                        inDic["CIS"]=lines[i][60:81]
                        i+=1
                        inDic["i0"]=lines[i][3:22]
                        inDic["Crc"]=lines[i][22:41]
                        inDic["omega"]=lines[i][41:60]
                        inDic["OMEGA DOT"]=lines[i][60:81]
                        i+=1
                        inDic["IDOT"]=lines[i][3:22]
                        inDic["L2"]=lines[i][22:41]
                        inDic["GPS Week"]=lines[i][41:60]
                        inDic["L2 P"]=lines[i][60:81]
                        i+=1
                        inDic["SV accuracy"]=lines[i][3:22]
                        inDic["SV health"]=lines[i][22:41]
                        inDic["TGD"]=lines[i][41:60]
                        inDic["IODC"]=lines[i][60:81]
                        i+=1
                        inDic["Transmission time"]=lines[i][3:22]
                        inDic["Fit interval"]=lines[i][22:41]

                          
                        
                        if name in dict_:
                            dict_[name][time]=inDic 
                        else:   
                            timeDic[time]=inDic
                            dict_[name]=timeDic

                    else:i+=1
        return dict_        
        
                    
                    
                    
                    
        #return dict object with key and dict of values
        
        
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
        
        
        
        
        
        
        
        
        