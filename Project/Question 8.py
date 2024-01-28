#Code MECH 309
#Question 8
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

def f_w1(x,q):
    f1 = np.sin(phi1(x)*q)
    return f1


N = int(input('Enter N = '))
L = 1
a = 0
b = 1
D = 1
q = np.arange(qstart,qend,0.1)
l = int((qend - qstart)/0.1 + 0.1)
f_n = 100
g = np.zeros(l)
f = np.zeros(l)
p_0 = qstart
p_1 = qend
list_of_p = np.zeros(l)
epsi = float(input('Enter espilon = '))
for i in range(l):
    g[i] = -(I_2_Trapezoid(a,b,N,q[i],L,f_n)/I_1_Trapezoid(a,b,N,q[i]))



def w(x):
    q = 6.8567
    h = float(b - a) / N
    intg = 0.0
    intg += f_w1(a,q)/2.0
    for i in range(int(1), int(N)):
        intg += f_w1(a + i*h,q)
    intg += f_w1(b,q)/2.0
    return intg * h


N = int(input('Enter N = '))
L = 1
a = 0
b = 1
D = 1
q = np.arange(qstart,qend,0.1)
l = int((qend - qstart)/0.1 + 0.1)
f_n = 100
g = np.zeros(l)
f = np.zeros(l)
p_0 = qstart
p_1 = qend
list_of_p = np.zeros(l)
epsi = float(input('Enter espilon = '))
for i in range(l):
    g[i] = -(I_2_Trapezoid(a,b,N,q[i],L,f_n)/I_1_Trapezoid(a,b,N,q[i]))
    while abs(p_1-p_0)> epsi:

    
