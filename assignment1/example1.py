

import numpy as np
import math

#Read File
file = open('GGM03S.txt','r') 
nlist = []
key = []
mlist = []
dict = {}

for line in file:
    line = line.strip()
    c = line[13:33]
    s = line[34:54]
    n = int(line[6:9])
    m = int(line[9:12])
    nlist += [n]
    mlist += [m]
    nm = str(n) + str(m)
    key += [nm]
    dict[int(nm)] = [float(c.replace('D','E')),float(s.replace('D','E'))]
file.close()

###Station 51 Test Data
##lat = np.radians(-24.72167)
##lon = np.radians(28.37667)

###Station 52 Test Data
##lat = np.radians(-24.68333)
##lon = np.radians(28.41667)


###Station 53 Test Data
##lat = np.radians(-24.66667)
##lon = np.radians(28.48500)

#Station 54 Test Data
lat = np.radians(-24.67500)
lon = np.radians(28.52500)


'_________________________________Constants__________________________________'

a = 6378137.0             
b = 6356752.3141         
GM = 3986005.0e8         
J2 = 108263.0e-8          
omega = 7292115.0e-11     
gammaE = 9.7803267715   
E = np.sqrt(a**2 - b**2)
e = E / a

g1 = 0.0052790414 *((np.sin(lat)))**2
g2 = 0.0000232718 *((np.sin(lat)))**4
g3 = 0.0000001262 *((np.sin(lat)))**6
g4 = 0.0000000007 *((np.sin(lat)))**8
gamma = gammaE * (1 + g1 + g2 + g3 + g4)    


R = a * (np.sqrt(1 - ((((e**2)*(1 - e**2))*((np.sin(lat))**2)) / \
                      ((1 - (e**2)*((np.sin(lat))**2))))))

P = {}
t = np.cos(((np.pi)/2)-(np.arctan(((b/a)**2)*np.tan(lat))))  
u = np.sin(((np.pi)/2)-(np.arctan(((b/a)**2)*np.tan(lat))))
P[00] = 1         
P[10] = t
P[11] = u
P[20] = (1.5)*(t**2) - (0.5)
P[21] = 3* u * t
P[22] = 3 * (u**2)

i = 86   #GGMOs file truncated at 86

'_________________________________Calculating P__________________________________'

for n in range(2,i):        
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
            h = t ** (n - m -2*k)
            
            this_loop = a* (c/(d*f*g))*h

            Sum += this_loop

        Pnm = (2 ** (-n) * (1 - t ** 2) ** (m/2.0)) * Sum
        nm = str(n) + str(m)
        P[int(nm)] = Pnm

J = {2: 108263.0e-8,4: -0.00000237091222,6: 0.00000000608347,8: -0.00000000001427}
Jn = {}

j = 2        
while j <= 8:
    n = j
    Jnorm = (np.sqrt(1/(2.0*n + 1.0))) * J[j] 
    Jn[j] = Jnorm
    j += 2

dict[20][0] = float(dict[20][0]) + Jn[2]
dict[40][0] = float(dict[40][0]) + Jn[4]
dict[60][0] = float(dict[60][0]) + Jn[6]
dict[80][0] = float(dict[80][0]) + Jn[8]

'_________________________________Normalize P__________________________________'

NormP = {}
for n in range(2,i):
    for m in range (0, n+1):
        nm = str(n) + str(m)
        pnm = P[int(nm)]

        if m == 0:
            j = 1
        if m!= 0:
            j = 2
        pn1 = j * (2.0*n+1.0)
        pn2 = float(math.factorial(n-m))
        pn3 = float(math.factorial(n+m))
        pn4 = pn2/pn3

        NormP[int(nm)] = (np.sqrt(pn1 * pn4)) * pnm
        
'_________________________________Calculate N__________________________________'

full_sum = 0
TotN = {}
for n in range (2, i):
    TotalN = 0
    a = 6378137.0 
    for m in range (0,n+1):
        
        nm = int(str(n) + str(m))
        C = dict[nm][0]
        S = dict[nm][1]
        TotalN += (C * np.cos(m*lon) + S * np.sin(m*lon)) * NormP[nm]
    TotalN = ((a/R)**n) * TotalN
    TotN[n] = TotalN
    
for n in range (2, i):
    full_sum = full_sum + TotN[n]
  

N = full_sum * (GM/(gamma * R))

print 'N 64 =', N       
            

    
    



























               
