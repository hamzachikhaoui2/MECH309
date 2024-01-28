#My Code
#MECH 309 Question 12
import matplotlib.pyplot as plt
import numpy as np
import random
from math import sin,cos,sinh,cosh,sqrt

########################## Functions from Q2.5.v2 ############################

def phi(x,r):
    phi = sin(r*x)+sinh(r*x)+((cos(r)+cosh(r))/(sin(r)+sinh(r)))*(cos(r*x)-cosh(r*x))
    return phi

def d2phi(x,r):
    d2phi=r**2*(sinh(r*x)-sin(r*x)+((cos(r)+cosh(r))/(sin(r)+sinh(r)))*(-cos(r*x)-cosh(r*x)))
    return d2phi

def f1(x,r,q1,q2):
    r1 = 1.87
    r2 = 4.69
    f1 = phi(x,r)*(d2phi(x,r1)*q1 + d2phi(x,r2)*q2)
    return f1

def f2(x,r,q1,q2):
    s = np.arange(x,1,0.01)
    n = int((1-x)/0.01)
    h = 0.01
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf+=(h/2)*(cos(q1*(phi(x,r1)-phi(s0,r1))+q2*(phi(x,r2)-phi(s0,r2)))+cos(q1*(phi(x,r1)-phi(s2,r1))+q2*(phi(x,r2)-phi(s2,r2))))
    f2 = 100*phi(x,r)*subf
    return f2

def numintg(a,b,n,x,q1,q2,f,r):   #numerical integration by trapezoidal rule with starting point (a), end point (b),
    intg = 0                #number of intervals (n), list of x values(x), variable q1 (q), function to integrate (f)
    h = (b-a)/n
    for i in range(n):
        x0 = x[i]
        x2 = x[i]+h
        intg += (h/2)*(f(x0,r,q1,q2)+f(x2,r,q1,q2)) #- (h**5/90)*1    #final coeff term ('1') replaces f''''(c)
    return intg


######################### Iterating code #######################

r1 = 1.87
r2 = 4.69

a = 0
b = 1
n = 10
h = (b-a)/n
x = np.arange(a,b,h)

print('Choose interval of interest for f(q1,q2) :')                 #lets the user choose the interval of interest for q1 and q2
qstart = int(input('Starting point = '))                          #starting point
qend = int(input('End point = '))                                 #end point
m = int(input('Number of intervals = ')) + 1                      #asks for the number of subintervals in the chosen interval

#q_1 = np.outer(np.linspace(qstart,qend,m), np.ones(m))              #prepares values of q1 as a horizontal array
#q_2 = np.outer(np.linspace(qstart,qend,m), np.ones(m))
#q_2 = q_1.copy().T                                                  #prepares values of q2 as a vertical array

#print(q_1)
#print(q_2)

#f_1 = np.zeros((m,m))
#f_2 = np.zeros((m,m))
#k = 0

#for i in range(m):
#    for j in range(m):
#        f_1[i,j] = numintg(a,b,n,x,q_1[i,j],q_2[i,j],f1,r1) + numintg(a,b,n,x,q_1[i,j],q_2[i,j],f2,r1)
#        f_2[i,j] = numintg(a,b,n,x,q_1[i,j],q_2[i,j],f1,r2) + numintg(a,b,n,x,q_1[i,j],q_2[i,j],f2,r2)
#    k += 1
    #print(k)



    
def indexfinder(x,y,q1,q2):
    xdif = 10
    ydif = 10
    index = 0
    xcoord = 0
    ycoord = 0
    for i in range(len(q1)):
        if abs(x-q1[i,i]) < xdif:
            xdif = abs(x-q1[i,i])
            xcoord = i
    for j in range(len(q1)):
        if abs(y-q2[j,j]) < ydif:
            ydif = abs(y-q2[j,j])
            ycoord = j

    index = [xcoord,ycoord]
    return index





###############  Broyden Method ################
tol = 0.001                                                   #sets the tolerance to three decimal places
compt = 0                                                          #sets the number of iterations to 0
norm = 1                                                       #sets the norm as a value greater than tol to activate while loop
x0y0 = [0,0]                                                   #sets the starting point to (0,0)
values = np.arange(-1,2,0.1)                                #prepares the set of values possible for the coordinates of the starting point
x0y0[0] = 2                               #sets the x coordinate as a random value between 0 and 2
x0y0[1] = -1  
                              #sets the y coordinate as a random value between 0 and 2
#start = x0y0.copy()
A0 = np.matrix(np.identity(2))

#print('N = '+str(n)+' : [x,y] = ['+str(x0y0[0])+','+str(x0y0[1])+']')   #prints the coordinates of the starting point
#plt.scatter(x0y0[0],x0y0[1])                                            #plots the starting point on a graph

convergence = True


while norm > tol:

    x1 = numintg(a,b,n,x,x0y0[0],x0y0[1],f1,r1) + numintg(a,b,n,x,x0y0[0],x0y0[1],f2,r1)
    y1 = numintg(a,b,n,x,x0y0[0],x0y0[1],f1,r2) + numintg(a,b,n,x,x0y0[0],x0y0[1],f2,r2)
    x1y1 = [x1,y1]
    y = np.dot(np.linalg.inv(-A0),(x1y1))
    x0y0 += y
    x1_2 = numintg(a,b,n,x,x0y0[0],x0y0[1],f1,r1) + numintg(a,b,n,x,x0y0[0],x0y0[1],f2,r1)
    y1_2 = numintg(a,b,n,x,x0y0[0],x0y0[1],f1,r2) + numintg(a,b,n,x,x0y0[0],x0y0[1],f2,r2)
    x1y1_2 = [x1_2,y1_2]
    delta_f = x1y1_2-x1y1

    
    A0 += ((delta_f-np.dot(A0,y))*y.T)/((x1y1[0] - x0y0[0])**2 + (x1y1[1] - x0y0[1])**2)
    #print(A0)
    compt += 1
    #print(n)
    norm = sqrt((x1 - x0y0[0])**2 + (y1 - x0y0[1])**2)

   #print(n)
    plt.scatter(x0y0[0],x0y0[1]) 
    plt.scatter(x1,y1)
    #print('N = '+str(n)+' : [x,y] = ['+str(x1)+','+str(y1)+'] ; norm = '+str(norm))
    print(x1y1)

    if compt == 100:
        convergence = False
        print('Broyden does not converge')
        break


#print(x1y1)
plt.plot()
plt.grid()
plt.show()
