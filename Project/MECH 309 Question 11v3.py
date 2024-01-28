#Code MECH 309
#Question 10
#My code Hamza Chikhaoui
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin,cos,sinh,cosh
import random

r1 = 1.87
r2 = 4.69


def phi1(x,r):
    phi1 = np.sin(r*x)+np.sinh(r*x)+((np.cos(r)+np.cosh(r))/(np.sin(r)+np.sinh(r)))*(np.cos(r*x)-np.cosh(r*x))
    return phi1

def d2phi1_1(x,r):
    d2phi1=(r**2)*(np.sinh(r*x)-np.sin(r*x)-((np.cos(r)+np.cosh(r))/(np.sin(r)+np.sinh(r)))*(np.cos(r*x)+np.cosh(r*x)))
    return d2phi1


def f1(x,q,r,r_d2,q_1,q_2):
    f1 = D*phi1(x,r)*(d2phi1_1(x,r_d2)*q)
    return f1

def f2(x,q_1,r1,q_2,r2,r):
    s = np.arange(x,1,0.01)
    n = int((1-x)/0.01)
    h = 0.01
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf += (h/2)*(cos(q_1*(phi1(x,r1)-phi1(s0,r1)))+cos(q_1*(phi1(x,r1)-phi1(s2,r1)))+(cos(q_2*(phi1(x,r2)-phi1(s0,r2)))+cos(q_2*(phi1(x,r2)-phi1(s2,r2)))))     #also try with (h/2) instead if it doesnt work
    f2 = 100*phi1(x,r)*subf  #try with a different variable that stays constant for subf if needed
    return f2

def numintg(a,b,n,x,q_1,q_2,f,r1,r2,r,r_d2):   #numerical integration by trapezoidal rule with starting point (a), end point (b),
    intg = 0                #number of intervals (n), list of x values(x), variable q1 (q), function to integrate (f)
    h = (b-a)/n
    for i in range(n):
        x0 = x[i]
        x2 = x[i]+h
        intg += (h/2)*(f(x0,q_1,r1,q_2,r2,r)+f(x2,q_1,r1,q_2,r2,r)) #- (h**5/90)*1    #final coeff term ('1') replaces f''''(c)
    return intg


a = 0
b = 1
n = 100
D = 1
L = 1
h = (b-a)/n
x = np.arange(a,b,h)
m = 101

print('Choose interval of interest for f(q1) :')
q_1_start = float(input('Starting point q1 = '))
q_1_end = float(input('End point q1 = '))
q_2_start = float(input('Starting point q2 = '))
q_2_end = float(input('End point q2 = '))

q_1 = np.outer(np.linspace(q_1_start,q_1_end,m), np.ones(m))
q_2 = np.outer(np.linspace(q_2_start,q_2_end,m), np.ones(m))
l_1 = int((q_1_end - q_1_start)/0.1 + 0.1)
l_2 = int((q_2_end - q_2_start)/0.1 + 0.1)


f_n = 100

g_1 = np.zeros((l_1,l_1))
g_2 = np.zeros((l_2,l_2))
#g = (g_1,g_2)
#q = (q_1,q_2)


#Def index string function
def indexfinder(r):
    dif = 0.5
    index = 0
    for i in range(len(q)):
        if abs(r-q[i]) < dif:
            dif = abs(r-q[i])
            index = i
    return index


#Add fixed point iteration on that

#Finding Values of q1 and q2 from functions g1 and g2
compteur = 0
for i in range(l_1):
    for j in range (l_2):
        g_1[i][j] = -(numintg(a,b,n,x,q_1[i][j],q_2[i][j],f1,r1,r2,r1,r2)+numintg(a,b,n,x,q_1[i][j],q_2[i][j],f2,r1,r2,r1,1))/(numintg(a,b,n,x,q_1[i][j],q_2[i][j],f1,r1,r2,r1,r1))
        g_2[i][j] = -(numintg(a,b,n,x,q_1[i][j],q_2[i][j],f1,r1,r2,r2,r1)+numintg(a,b,n,x,q_1[i][j],q_2[i][j],f2,r1,r2,r2,1))/(numintg(a,b,n,x,q_1[i][j],q_2[i][j],f1,r1,r2,r2,r2))
        q_1[i][j] = g_1[i][j]
        q_2[i][j] = g_2[i][j]
    compteur += 1
    print(compteur)


val_1 = []
val_2 = []

for k in range(l_1):
    for l in range (l_2):
        val_1.append(g_1[k][l]-q_1[k][l])
        val_2.append(g_2[k][l]-q_2[k][l])

q_1_f = np.min(val_1)
q_2_f = np.min(val_2)
print(q_1_f,q_2_f)
