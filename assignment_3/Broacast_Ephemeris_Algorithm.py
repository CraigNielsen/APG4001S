

__author__ = 'Craig FERGUSON'
from numpy import *

class broadCastEphemeris(object):
    '''
    object for Broadcast Ephemeris Calculations : GEODESY
    '''

    def __init__(self, rinex_):
        '''
        Constructor
        '''
        self.rinex = rinex_
        self.A = None
        self.n0 = None
        self.tk = None
        self.n = None
        self.mk = None
        self.vk = None
        self.Ek = None
        self.latK = None
        self.duk = None
        self.drk = None
        self.dik = None
        self.uk = None
        self.rk = None
        self.ik = None
        self.xk_ = None
        self.yk_ = None
        self.omegK = None
        self.xk = None
        self.yk = None
        self.zk = None
        self.mu = 3.986008e14

    def setA(self, name,time):
        self.A = (self.rinex[name][time]["sqrtA"])**2

    def _n0(self):
        n0 = sqrt(self.mu / (self.A ** 3))
        self.n0 = n0
        return n0

    def _tk(self,name,time,t):
        tk = t - self.rinex[name][time]["Toe"]
        if tk > 302400: tk = tk - 604800
        if tk < -302400: tk = tk + 604800
        self.tk=tk
        return tk

    def _n(self,name,time):
        # corrected mean motion
        n = self.n0 + self.rinex[name][time]["DeltaN"]
        self.n=n
        return n

    def _mk(self,name,time):
        mk = self.rinex[name][time]["MO"] + self.n *self.tk
        # mk = self.rinex[name][time]["MO"] - self.n * self.tk
        self.mk = mk
        return mk

    def _Ek(self,name,time):
        # Keplers equation for anomaly
        Ek = self.mk + self.rinex[name][time]["e"] * sin(self.mk)
        i = 0
        while i < 4:
            Ek = self.mk + self.rinex[name][time]["e"] * sin(Ek)
            i += 1
        self.Ek = Ek
        return Ek

    def _vk(self,name,time):
        # true anomoly
        top = (sqrt(1 - (self.rinex[name][time]["e"] ** 2)) * sin(self.Ek)) / (1 - self.rinex[name][time]["e"] * cos(self.Ek))
        bottom = (cos(self.Ek) - self.rinex[name][time]["e"]) / (1 - self.rinex[name][time]["e"] * cos(self.Ek))
        vk = arctan(top / bottom)
        self.vk = vk
        return vk

    def _latk(self,name,time):
        # argument of latitude
        latk = self.vk + self.rinex[name][time]["omega"]
        self.latK = latk
        return latk

    def _duk(self,name,time):
        Cus =self.rinex[name][time]["Cus"]
        Cuc =self.rinex[name][time]["Cuc"]
        Ok = self.latK
        self.duk = Cus * sin(2 * Ok) + Cuc * cos(2 * Ok)
        return self.duk

    def _drk(self,name,time):
        Crc =   self.rinex[name][time]["Crc"]
        Crs =self.rinex[name][time]["Crs"]
        Ok = self.latK
        self.drk = Crc * cos(2 * Ok) + Crs * sin(2 * Ok)
        return self.drk

    def _dik(self,name,time):
        Cic =self.rinex[name][time]["Cic"]
        Cis =self.rinex[name][time]["CIS"]
        Ok = self.latK
        self.dik = Cic * cos(2 * Ok) + Cis * sin(2 * Ok)
        return self.dik

    def _uk(self):
        self.uk = self.latK + self.duk
        return self.uk

    def _rk(self,name,time):
        self.rk = self.A * (1 - self.rinex[name][time]["e"] * cos(self.Ek)) + self.drk
        return self.rk

    def _ik(self,name,time):
        i0 =self.rinex[name][time]["i0"]
        IDOT =self.rinex[name][time]["IDOT"]
        self.ik = i0 + self.dik + (IDOT) * self.tk
        return self.ik

    def _xk_(self):
        self.xk_ = self.rk * cos(self.uk)
        return self.xk_

    def _yk_(self):
        self.yk_ = self.rk * sin(self.uk)
        return self.yk_

    def _omegK(self,name,time):
        sigDOTe = 7.292115167e-5
        omega0=self.rinex[name][time]["OMEGA"]
        omegaDOT=self.rinex[name][time]["OMEGA DOT"]
        self.omegK = omega0 + (omegaDOT - sigDOTe) * self.tk - sigDOTe * self.rinex[name][time]["Toe"]
        return self.omegK

    def getxk(self):
        ok = self.omegK
        xk = self.xk_ * cos(ok) - self.yk_ * cos(self.ik) * sin(ok)
        return xk

    def getyk(self):
        ok = self.omegK
        yk = self.xk_ * sin(ok) + self.yk_ * cos(self.ik) * cos(ok)
        return yk

    def getzk(self):
        zk = self.yk_ * sin(self.ik)
        return zk

    def getEphemeris(self,name,time):
        if self.rinex == None:
            print "no rinex file read yet, please read rinex file"
            return
        ephemeris=[]
        ran=[]
        for t in range(0,96):
            t=t*900
            self.setA(name,time)
            self._n0()
            self._tk(name,time,t)
            self._n(name,time)
            self._mk(name,time)
            self._Ek(name,time)
            self._vk(name,time)
            self._latk(name,time)
            self._duk(name,time)
            self._drk(name,time)
            self._dik(name,time)
            self._uk()
            self._rk(name,time)
            self._ik(name,time)
            self._xk_()
            self._yk_()
            self._omegK(name,time)
            x= self.getxk()
            y= self.getyk()
            z= self.getzk()
            ephemeris.append([x,y,z])
            r=sqrt(x**2 + y**2 + z**2)
            ran.append(r)
        return ephemeris,ran