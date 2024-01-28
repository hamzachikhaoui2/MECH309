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
q1_1 = 4.24
q1_2 = 1.72
q2_2 = -0.83
q1_3 = 1.775
q2_3 = -0.833
q3_3 = -0.0278

def psi_1(x):
    return q1_1*phi(x,r1)
def psi_2(x):
    return q1_2*phi(x,r1)+q2_2*phi(x,r2)
def psi_3(x):
    return q1_3*phi(x,r1)+q2_3*phi(x,r2)+q3_3*phi(x,r3)

a = np.arange(0,1,0.00001)

y_1 = [0]*len(a)
y_2 = [0]*len(a)
y_3 = [0]*len(a)


for i in range(len(a)):
    y_1[i] = psi_1(a[i])
    y_2[i] = psi_2(a[i])
    y_3[i] = psi_3(a[i])

plt.plot(a,y_1,'r-',label = 'psi1 with q1')
plt.plot(a,y_2,'b-',label = 'psi2 with q1 and q2')
plt.plot(a,y_3,'g-',label = 'psi3 with q1 and q2 and q3')
plt.grid()
plt.title("Convergence of psi_1,psi_2,psi_3")
plt.legend()
plt.show()






