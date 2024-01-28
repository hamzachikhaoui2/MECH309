import matplotlib.pyplot as plt
import numpy as np
import random
from math import sin,cos,sinh,cosh,sqrt

########################## Functions ############################

def phi(x,r):
    phi = sin(r*x)+sinh(r*x)+((cos(r)+cosh(r))/(sin(r)+sinh(r)))*(cos(r*x)-cosh(r*x))
    return phi

def d2phi(x,r):
    d2phi=r**2*(sinh(r*x)-sin(r*x)+((cos(r)+cosh(r))/(sin(r)+sinh(r)))*(-cos(r*x)-cosh(r*x)))
    return d2phi

def f1(x,r,q1,q2,q3):
    r1 = 1.87
    r2 = 4.69
    r3 = 7.85
    f1 = phi(x,r)*(d2phi(x,r1)*q1 + d2phi(x,r2)*q2 + d2phi(x,r3)*q3)
    return f1

def f2(x,r,q1,q2,q3):
    s = np.arange(x,1,0.01)
    n = int((1-x)/0.01)
    h = 0.01
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf+=(h/2)*(cos(q1*(phi(x,r1)-phi(s0,r1))+q2*(phi(x,r2)-phi(s0,r2))+q3*(phi(x,r3)-phi(s0,r3)))+cos(q1*(phi(x,r1)-phi(s2,r1))+q2*(phi(x,r2)-phi(s2,r2))+q3*(phi(x,r3)-phi(s2,r3))))
    f2 = 100*phi(x,r)*subf
    return f2

def numintg(a,b,n,x,q1,q2,q3,f,r):   #numerical integration by trapezoidal rule with starting point (a), end point (b),
    intg = 0                #number of intervals (n), list of x values(x), variable q1 (q), function to integrate (f)
    h = (b-a)/n
    for i in range(n):
        x0 = x[i]
        x2 = x[i]+h
        intg += (h/2)*(f(x0,r,q1,q2,q3)+f(x2,r,q1,q2,q3)) #- (h**5/90)*1    #final coeff term ('1') replaces f''''(c)
    return intg

######################### Iterating code #######################

r1 = 1.87
r2 = 4.69
r3 = 7.85

a = 0
b = 1
n = 10
h = (b-a)/n
x = np.arange(a,b,h)

tol = float(input('Choose tolerance : '))

norm = 1
x0y0z0 = [0,0,0]                                                   #sets the starting point to (0,0)
values = np.arange(-1,2,0.01)                                  #prepares the set of values possible for the coordinates of the starting point
x0y0z0[0] = random.choice(values)                                #sets the x coordinate as a random value between 0 and 2
x0y0z0[1] = random.choice(values)
x0y0z0[2] = random.choice(values)                                #sets the y coordinate as a random value between 0 and 2
A0 = np.matrix(np.identity(3))

k = 0

while norm > tol:
    f1x0 = numintg(a,b,n,x,x0y0z0[0],x0y0z0[1],x0y0z0[2],f1,r1) + numintg(a,b,n,x,x0y0z0[0],x0y0z0[1],x0y0z0[2],f2,r1)
    f2x0 = numintg(a,b,n,x,x0y0z0[0],x0y0z0[1],x0y0z0[2],f1,r2) + numintg(a,b,n,x,x0y0z0[0],x0y0z0[1],x0y0z0[2],f2,r2)
    f3x0 = numintg(a,b,n,x,x0y0z0[0],x0y0z0[1],x0y0z0[2],f1,r3) + numintg(a,b,n,x,x0y0z0[0],x0y0z0[1],x0y0z0[2],f2,r3)
    fx0 = np.matrix([[f1x0],[f2x0],[f3x0]])
    
    yvector = (A0.I)*(-1*fx0)
    
    f1x0plusy = numintg(a,b,n,x,x0y0z0[0]+yvector[0,0],x0y0z0[1]+yvector[1,0],x0y0z0[2]+yvector[2,0],f1,r1) + numintg(a,b,n,x,x0y0z0[0]+yvector[0,0],x0y0z0[1]+yvector[1,0],x0y0z0[2]+yvector[2,0],f2,r1)
    f2x0plusy = numintg(a,b,n,x,x0y0z0[0]+yvector[0,0],x0y0z0[1]+yvector[1,0],x0y0z0[2]+yvector[2,0],f1,r2) + numintg(a,b,n,x,x0y0z0[0]+yvector[0,0],x0y0z0[1]+yvector[1,0],x0y0z0[2]+yvector[2,0],f2,r2)
    f3x0plusy = numintg(a,b,n,x,x0y0z0[0]+yvector[0,0],x0y0z0[1]+yvector[1,0],x0y0z0[2]+yvector[2,0],f1,r3) + numintg(a,b,n,x,x0y0z0[0]+yvector[0,0],x0y0z0[1]+yvector[1,0],x0y0z0[2]+yvector[2,0],f2,r3)
    fx0plusy = np.matrix([[f1x0plusy],[f2x0plusy],[f3x0plusy]])
    
    deltaf = fx0plusy - fx0
    x0y0z0 = [x0y0z0[0]+yvector[0,0],x0y0z0[1]+yvector[1,0],x0y0z0[2]+yvector[2,0]]
    ynorm = sqrt(yvector[0,0]**2 + yvector[1,0]**2+ yvector[2,0]**2)
    A0 = A0 + ((deltaf - A0*yvector)*(yvector.T))/(ynorm**2)
    norm = ynorm
    k += 1
    #print('n = '+str(k)+', yvector = '+str(yvector))
    #print('x0 = '+str(x0y0))

print('q1 = '+str(x0y0z0[0])+', q2 = '+str(x0y0z0[1]),', q3 = '+str(x0y0z0[2]))
