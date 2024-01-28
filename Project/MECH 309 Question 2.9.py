#!/usr/bin/env python
# coding: utf-8

# In[26]:


import matplotlib.pyplot as plt
import numpy as np
import random
from math import sin,cos,sinh,cosh


# In[39]:


########################## Function #############################

def phi1(x):
    phi1 = sin(1.87*x)+sinh(1.87*x)+((cos(1.87)+cosh(1.87))/(sin(1.87)+sinh(1.87)))*(cos(1.87*x)-cosh(1.87*x))
    return phi1

def d2phi1(x):
    d2phi1=1.87**2*(sinh(1.87*x)-sin(1.87*x)+((cos(1.87)+cosh(1.87))/(sin(1.87)+sinh(1.87)))*(-cos(1.87*x)-cosh(1.87*x)))
    return d2phi1

def f1(x,q,f_n):
    f1 = phi1(x)*d2phi1(x)
    return f1

def f2(x,q,f_n):
    s = np.arange(x,1,0.1)
    n = int((1-x)/0.1)
    h = 0.1
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf += (h/2)*(cos(q*(phi1(x)-phi1(s0)))+cos(q*(phi1(x)-phi1(s2))))     #also try with (h/2) instead if it doesnt work
    f2 = phi1(x)*subf*f_n                                                       #try with a different variable that stays constant for subf if needed
    return f2

def numintg(a,b,n,x,q,f,f_n):   #numerical integration by trapezoidal rule with starting point (a), end point (b),
    intg = 0                #number of intervals (n), list of x values(x), variable q1 (q), function to integrate (f)
    h = (b-a)/n
    for i in range(n):
        x0 = x[i]
        x2 = x[i]+h
        intg += (h/2)*(f(x0,q,f_n)+f(x2,q,f_n)) #- (h**5/90)*1   #final coeff term ('1') replaces f''''(c)
    return intg


# In[40]:


######################### Iterating code #######################

a = 0
b = 1
n = 100
h = (b-a)/n
x = np.arange(a,b,h)
L = 1

print('Choose interval of interest for f(q1) :')
qstart = float(input('Starting point = '))
qend = float(input('End point = '))
q = np.arange(qstart,qend,h)
phi1_L = phi1(L)

f_n = np.arange(0,101,1)




l = int((qend - qstart)/h + h)
psih_L = []

f = np.zeros(l)





def indexfinder(r,q):
    diff = np.abs(q - r)
    index = np.where(diff == min(diff))
    return q[index]


# In[50]:


############################## Secant method with varying f_n ################################

for i in range (len(f_n)):
    for j in range(len(f_n)):
        f[j] = ((q[j])*(numintg(a,b,n,x,q[j],f1,1))) + (numintg(a,b,n,x,q[j],f2,f_n[i]))
    def indexfinder(r):
        dif = 0.5
        index = 0
        for i in range(len(q)):
            if abs(r-q[i]) < dif:
                dif = abs(r-q[i])
                index = i
        return index

    r0 = (random.randint(400,450))/100
    r1 = (random.randint(400,450))/100
    while r1 == r0:
        r1 = (random.randint(400,450))/100
    f0 = f[indexfinder(r0)]
    f1v = f[indexfinder(r1)]
    if float(f1v-f0) != 0 :
        ri = r1 - f1v*((r1-r0)/(float(f1v-f0)))
    else:
        print("f1v and f0 are the same")
    tol = 0.005                                      #lets the user choose the tolerance
    while abs(ri-r1) >= tol:
        r0list = [r1]
        r1list = [ri]
        r0 = r0list[0]
        r1 = r1list[0]
        f0list = [f1v]
        f0 = f0list[0]
        f1v = f[indexfinder(r1)]
        ri = r1 - f1v * ((r1-r0)/float((f1v-f0)))

    q1str = str((ri)*100).split('.',1)
    q1 = float(q1str[0])/100

    psih_L.append(q1*phi1_L)


# In[53]:


len(psih_L)


# In[66]:


fig,ax = plt.subplots(1,1,figsize=(10,6))

ax.set_xlabel('f_n')
ax.set_ylabel('psih_L')
ax.plot(f_n, psih_L[:101], color = 'blue',label = "psih_L")
ax.legend()
ax.set_title('Plot for Question 2.9')

plt.grid()
plt.show()


# In[ ]:





