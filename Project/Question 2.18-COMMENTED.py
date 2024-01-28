import matplotlib.pyplot as plt
import numpy as np
import random
from math import sin,cos,sinh,cosh,sqrt
#from scipy.signal import savgol_filter

########################## Functions ############################

def phi(x,r):
    phi = sin(r*x)+sinh(r*x)+((cos(r)+cosh(r))/(sin(r)+sinh(r)))*(cos(r*x)-cosh(r*x)) #Function of Phi(x) from (4) in the report
    return phi

def d2phi(x,r):
    d2phi=r**2*(sinh(r*x)-sin(r*x)+((cos(r)+cosh(r))/(sin(r)+sinh(r)))*(-cos(r*x)-cosh(r*x))) #Function of the double derivative of the 
                                                                                                #function of Phi(x) from (4) in the report
    return d2phi

def f1(x,r,q1,q2,q3): #Function f1 that was used in previous integration questions. Not used for the iteration loop of this code.
    r1 = 1.87
    r2 = 4.69
    r3 = 7.85
    f1 = phi(x,r)*(d2phi(x,r1)*q1 + d2phi(x,r2)*q2 + d2phi(x,r3)*q3)
    return f1



def f2_1(x,q1,q2,q3): #Integration function from x to L = 1 for Psi1 = q1*phi1
    s = np.arange(x,1,0.01)
    n = int((1-x)/0.01)
    h = 0.01
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf += (h/2)*(cos(q1*(phi(x,r1)-phi(s0,r1)))+cos(q1*(phi(x,r1)-phi(s2,r1))))     #Using Trapezoidal Rule
                                                                         #also try with (h/2) instead if it doesnt work
    f21 = 100*subf    #FN incorporated into integration function                     #try with a different variable that stays constant for subf if needed
    return f21



def f2_2(x,q1,q2,q3):     #Integration function from x to L = 1 for Psi2 = q1*phi1+q2*phi2
    s = np.arange(x,1,0.01)
    n = int((1-x)/0.01)
    h = 0.01
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf+= (h/2)*(cos(q1*(phi(x,r1)-phi(s0,r1))+q2*(phi(x,r2)-phi(s0,r2)))  #Using Trapezoidal Rule
                      +cos(q1*(phi(x,r1)-phi(s2,r1))+q2*(phi(x,r2)-phi(s2,r2))))
    f22 = 100*subf  #FN incorporated into integration function
    return f22

def f2_3(x,q1,q2,q3):   #Integration function from x to L = 1 for Psi3 = q1*phi1+q2*phi2 +q3*phi3
    s = np.arange(x,1,0.01)
    n = int((1-x)/0.01)
    h = 0.01
    subf = 0
    for i in range(n):
        s0 = s[i]
        s2 = s[i]+h
        subf+=(h/2)*(cos(q1*(phi(x,r1)-phi(s0,r1))+q2*(phi(x,r2)-phi(s0,r2))+q3*(phi(x,r3)-phi(s0,r3))) #Using Trapezoidal Rule
                     +cos(q1*(phi(x,r1)-phi(s2,r1))+q2*(phi(x,r2)-phi(s2,r2))+q3*(phi(x,r3)-phi(s2,r3))))
    f23 = 100*subf  #FN incorporated into integration function
    return f23

def numintg(a,b,n,x,q1,q2,q3,f):   #numerical integration by trapezoidal rule with starting point (a), end point (b),
    intg = 0                #number of intervals (n), list of x values(x), variable q1 (q), function to integrate (f)
    h = (b-a)/n
    for j in range(n):
        x0 = x
        x2 = x+h
        intg += (h/2)*(f(x0,q1,q2,q3)+f(x2,q1,q2,q3)) #- (h**5/90)*1    #final coeff term ('1') replaces f''''(c)
    return intg



######################### Useful function 2 #######################

r1 = 1.87 #Root 1
r2 = 4.69  #Root 2
r3 = 7.85  #Root 3
q1_1 = 4.24 #q1 obtained with Psi1
q1_2 = 1.72 # q1 obtained with Psi2
q2_2 = -0.83 #q2 obtained with Psi2
q1_3 = 1.775 # q1 obtained for Psi3
q2_3 = -0.833  #q2 obtained for Psi3
q3_3 = -0.0278 # q3 obtained for Psi3

def psi_1(x):        #Psi1 Definition
    return q1_1*phi(x,r1)
def psi_2(x):
    return q1_2*phi(x,r1)+q2_2*phi(x,r2) #Psi2 Definition
def psi_3(x):
    return q1_3*phi(x,r1)+q2_3*phi(x,r2)+q3_3*phi(x,r3) #Psi3 Definition

def d_psi_1(x):   #Double derivative with respect to x of Psi 1
    return q1_1*d2phi(x,r1)
def d_psi_2(x):   #Double derivative with respect to x of Psi 2
    return q1_2*d2phi(x,r1)+q2_2*d2phi(x,r2) 
def d_psi_3(x):   #Double derivative with respect to x of Psi 3
    return q1_3*d2phi(x,r1)+q2_3*d2phi(x,r2)+q3_3*d2phi(x,r3)


def f_inter_plot(f,x): #Interpolation function with degree of integrating polynomial = 5.
    return f[5]+f[4]*x+f[3]*x**2+f[2]*x**3+f[1]*x**4+f[0]*x**5



############################## Iterating code ###############################


a = np.arange(0,1.01,0.01)
b = np.arange(0,1.01,0.001)
D = 1 #D = 1 added automatically to integrating functions
FN = 100 #Value of FN already added previously to integrating functions

f_1 = np.zeros(101) #solution of the governing equation (2) using Psi1
f_2 = np.zeros(101) #solution of the governing equation (2) using Psi2
f_3 = np.zeros(101) #solution of the governing equation (2) using Psi3


for i in range(len(a)): #Computing the solution of the governing equation (2) using Psi1, Psi2, and Psi3.
    f_1[i] = (d_psi_1(a[i]))+ numintg(0,1,100,a[i],q1_1,1,1,f2_1)
    f_2[i] = (d_psi_2(a[i]))+ numintg(0,1,100,a[i],q1_2,q2_2,1,f2_2)
    f_3[i] = (d_psi_3(a[i]))+ numintg(0,1,100,a[i],q1_3,q2_3,q3_3,f2_3)





plt.plot(a,f_1,'r-',label = 'Equation (2) with psi1 with q1') #Plotting the obtained solutions of the governing equation (2) using Psi1, Psi2, and Psi3.
plt.plot(a,f_2,'b-',label = 'Equation (2) with psi2 with q1 and q2')
plt.plot(a,f_3,'g-',label = 'Equation (2) with psi3 with q1 and q2 and q3')


plt.grid()
plt.title("Convergence of psi_1,psi_2,psi_3 [1]")
plt.legend()
plt.show()





########################################## Plot of psi1 psi2 psi3 ###############################################

y_1 = [0]*len(a)
y_2 = [0]*len(a)
y_3 = [0]*len(a)


for i in range(len(a)): #Plotting the functions Psi1, Psi2, and Psi3 to check for some convergence in their behavior
    y_1[i] = psi_1(a[i])
    y_2[i] = psi_2(a[i])
    y_3[i] = psi_3(a[i])

plt.plot(a,y_1,'r-',label = 'psi1 with q1')
plt.plot(a,y_2,'b-',label = 'psi2 with q1 and q2')
plt.plot(a,y_3,'g-',label = 'psi3 with q1 and q2 and q3')


plt.grid()
plt.title("Convergence of psi_1,psi_2,psi_3 [2] ")
plt.legend()
plt.show()



###################################### Polynomial Interpolation Regular #####################################################



f_1_polyf = np.polyfit(a,f_1,5) #Obtain the polynomial interpolations of the solutions of the governing equation (2) using a degree 5.
f_2_polyf = np.polyfit(a,f_2,5)
f_3_polyf = np.polyfit(a,f_3,5)



f_1_inter_pl = np.zeros(len(a))
f_2_inter_pl = np.zeros(len(a))
f_3_inter_pl = np.zeros(len(a))




for i in range (len(a)):
    f_1_inter_pl[i] = f_inter_plot(f_1_polyf,a[i]) #Put a sequence of points through the interpolation equation, and save the obtained points
    f_2_inter_pl[i] = f_inter_plot(f_2_polyf,a[i])
    f_3_inter_pl[i] = f_inter_plot(f_3_polyf,a[i])
    

plt.plot(a,f_1_inter_pl,'r-',label = 'polynomial fit of Equation (2) profile with psi1 with q1') #Plot the obtained points to model the polynomial interpolation 
plt.plot(a,f_2_inter_pl,'b-',label = 'polynomial fit of Equation (2) profile with psi2 with q1 and q2')
plt.plot(a,f_3_inter_pl,'g-',label = 'polynomial fit of Equation (2) profile with psi3 with q1 and q2 and q3')

plt.grid()
plt.title("Convergence of psi_1,psi_2,psi_3 [3] ")
plt.legend()
plt.show()





#################################################### Polynomial Interpolation with filter ######################################################################

f_1_inter = savgol_filter(f_1, 51, 3) #Use of a Savgol filter to confirm the results obtained from the polynomial interpolation.
f_2_inter = savgol_filter(f_2, 51, 3)
f_3_inter = savgol_filter(f_3, 51, 3)

plt.plot(a,f_1_inter,'r-',label = 'Filtered psi1 with q1')
plt.plot(a,f_2_inter,'b-',label = 'Filtered psi2 with q1 and q2')
plt.plot(a,f_3_inter,'g-',label = 'Filtered psi3 with q1 and q2 and q3')

plt.grid()
plt.title("Convergence of psi_1,psi_2,psi_3 [4] ")
plt.legend()
plt.show()



