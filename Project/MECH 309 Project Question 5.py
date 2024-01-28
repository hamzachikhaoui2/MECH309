#Code MECH 309
#Question 6
#My code Hamza Chikhaoui
import numpy as np
import matplotlib.pyplot as plt
import math

r0 = 1.875


def phi1(x):
    phi1 = np.sin(r0*x)+np.sinh(r0*x)+((np.cos(r0)+np.cosh(r0))/(np.sin(r0)+np.sinh(r0)))*(np.cos(r0*x)-np.cosh(r0*x))
    return phi1

def d2phi1_1(x):
    d2phi1=(r0**2)*(np.sinh(r0*x)-np.sin(r0*x)-((np.cos(r0)+np.cosh(r0))/(np.sin(r0)+np.sinh(r0)))*(np.cos(r0*x)+np.cosh(r0*x)))
    return d2phi1


def f1(x,q,D):
    f1 = phi1(x)*d2phi1_1(x)*q*D
    return f1


def I_1_Trapezoid(a,b,N,q):

    h = float(b - a) / N
    intg = 0.0
    intg += f1(a,q,D)/2.0
    for i in range(int(1), int(N)):
        intg += f1(a + i*h,q,D)
    intg += f1(b,q,D)/2.0
    return intg * h
    
    

def f2(x,q,L): 
    f2 = ((phi1(x)*np.cos(float(q)*(phi1(x)-phi1(L))))+phi1(x))*float(((L-x)/2))
    return f2


#def f2(x,q,L):
    s = np.arange(x,1,0.1)
    n = int((1-x)/0.1)
    h = 0.1
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf += h*(np.cos(float(q)*(phi1(x)-phi1(s0)))+np.cos(float(q)*(phi1(x)-phi1(s2))))
    f2 = phi1(x)*subf
    return f2

    
#def f2_inter_Trapez(L,x,x1,x2,q):
#    f2_inter = ((x1-x2)/2)*(f2(x,x,x1,x2,q)+f2(x,L,x1,x2,q)) #Only s changes in here
#    return f2_inter


def I_2_Trapezoid(a,b,N,q,L,f_n):
    h = float(b - a) / N
    intg_2 = 0.0
    intg_2 += f2(a,q,L)/2.0
    for i in range(1, N):
        intg_2 += f2(a + i*h,q,L)
    intg_2 += f2(b,q,L)/2.0
    return intg_2 * h * f_n

print('Choose interval of interest for f(q1) :')
qstart = float(input('Starting point = '))
qend = float(input('End point = '))

N = int(input('Enter N = '))
L = 1
a = 0
b = 1
D = 1
q = np.arange(qstart,qend,0.1)
l = int((qend - qstart)/0.1 + 0.1)
f_n = 100
f = np.zeros(l)

for i in range(l):
    f[i] = I_1_Trapezoid(a,b,N,q[i]) + I_2_Trapezoid(a,b,N,q[i],L,f_n)
 


#ftest = numintg(a,b,n,x,1,f1)3
#print(ftest)

plt.plot(q,f)
plt.title('f(q1)')
plt.xlabel('q1')
plt.ylabel('f')
plt.grid()
plt.show()

    
