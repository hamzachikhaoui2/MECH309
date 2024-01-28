#Code MECH 309
#Question 10
#My code Hamza Chikhaoui
import numpy as np
import matplotlib.pyplot as plt
import math

r1 = 1.87
r2 = 4.69


def phi1(x,r):
    phi1 = np.sin(r*x)+np.sinh(r*x)+((np.cos(r)+np.cosh(r))/(np.sin(r)+np.sinh(r)))*(np.cos(r*x)-np.cosh(r*x))
    return phi1

def d2phi1_1(x,r):
    d2phi1=(r**2)*(np.sinh(r*x)-np.sin(r*x)-((np.cos(r)+np.cosh(r))/(np.sin(r)+np.sinh(r)))*(np.cos(r*x)+np.cosh(r*x)))
    return d2phi1


def f1(x,D,q_1,r1,q_2,r2,r):
    f1 = D*phi1(x,r)*((d2phi1_1(x,r1)*q_1)+(d2phi1_1(x,r2)*q_2))
    return f1

def f2(x,q):
    s = np.arange(x,1,0.01)
    n = int((1-x)/0.01)
    h = 0.01
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf += (h/2)*(cos(q*(phi1(x)-phi1(s0)))+cos(q*(phi1(x)-phi1(s2))))     #also try with (h/2) instead if it doesnt work
    f2 = 100*phi1(x)*subf                                                       #try with a different variable that stays constant for subf if needed
    return f2

def numintg(a,b,n,x,q,f):   #numerical integration by trapezoidal rule with starting point (a), end point (b),
    intg = 0                #number of intervals (n), list of x values(x), variable q1 (q), function to integrate (f)
    h = (b-a)/n
    for i in range(n):
        x0 = x[i]
        x2 = x[i]+h
        intg += (h/2)*(f(x0,q)+f(x2,q)) #- (h**5/90)*1    #final coeff term ('1') replaces f''''(c)
    return intg


a = 0
b = 1
n = 100
D = 1
L = 1
h = (b-a)/n
x = np.arange(a,b,h)


print('Choose interval of interest for f(q1) :')
q_1_start = float(input('Starting point q1 = '))
q_1_end = float(input('End point q1 = '))
q_2_start = float(input('Starting point q2 = '))
q_2_end = float(input('End point q2 = '))

q_1 = np.arange(q_1_start,q_1_end,0.1)
q_2 = np.arange(q_2_start,q_2_end,0.1)

l_1 = int((q_1_end - q_1_start)/0.1 + 0.1)
l_2 = int((q_2_end - q_2_start)/0.1 + 0.1)


f_n = 100
f_1 = np.zeros(l_1)
f_2 = np.zeros(l_2)



### The plot part starts here. For some reason, it wouldnt plot.
### Plot the thing in 3d? (with i's and j's??)
### Plot the thing in 2D, simply with q1 and q2 as changing variables "curseurs"??)
for i in range(l_1):
    for j in range (l_2):
        f_1[i] = I_1_Trapezoid(a,b,N,q_1[i],r1,q_2[j],r2,D,r1) + I_2_Trapezoid(a,b,N,q_1[i],q_2[j],L,f_n,r1,r1,r2)
        f_2[i] = I_1_Trapezoid(a,b,N,q_1[i],r2,q_2[j],r2,D,r2) + I_2_Trapezoid(a,b,N,q_1[i],q_2[j],L,f_n,r2,r1,r2)
    


#ftest = numintg(a,b,n,x,1,f1)3
#print(ftest)

    

#q_1,q_2 = np.meshgrid(q_1,q_2)
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#ax.plot_surface(q_1, q_2, f_1,rstride=8, cstride=8, alpha=0.2)
#ax.plot_surface(q_1, q_2, f_2,rstride=8, cstride=8, alpha=0.2)
#cset = ax.contour(q_1, q_2, f_1, zdir='z', offset=0, cmap=cm.coolwarm)
#cset = ax.contour(q_1, q_2, f_2, zdir='z', offset=0, cmap=cm.coolwarm)


plt.plot(q_1,f_1)
#plt.plot(q_2,f_2)
plt.title('f(qi)')
plt.xlabel('qi')
plt.ylabel('fi')
plt.legend()
plt.grid(True)
plt.show()
