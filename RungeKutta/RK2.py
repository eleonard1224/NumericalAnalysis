'''
Implements 2nd-order Runge-Kutta (RK2) method for solving a first-order ordinary diff eqn
'''

# Importing the libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Defining the RK2 method
# y'(t) = f(t,y(t))
# y0 = y(0)
# dt is a fixed time step size
# ti is initial time
# tf is final time
# RK2 returns an array filled with values of y
#     and an array filled with values of t
def RK2(f,y0,dt,ti,tf):
    
    tvals = np.array([]);
    yvals = np.array([]);
    
    t = ti;
    y = y0;
    tvals = np.append(tvals, t);
    yvals = np.append(yvals, y);
  
    while t < tf:
        t = t + dt;
        k1 = dt*f(y,t);
        k2 = dt*f(y+k1,t+dt);
        y = y + 0.5*(k1+k2);
        tvals = np.append(tvals, t);
        yvals = np.append(yvals, y);
    
    return tvals, yvals;


# Testing the solution to y'(t) = -y(t), y(0) = 1 => y(t) = e^(-t)
y0 = 1.0; dt = 0.01; ti = 0.0; tf = 1.0; 
def f1(y,t):
    return -1.0*y;
t1, y1 = RK2(f1,y0,dt,ti,tf);
tref1 = np.linspace(ti,tf,100);
yref1 = np.exp(-1.0*tref1);

# Plotting
mpl.rcParams.update({'font.size': 16})
# mpl.rc('lines', linewidth=4)

fig1 = plt.figure(1)
plt.plot(tref1,yref1,'r',label='ref')
plt.plot(t1[0:-1:10],y1[0:-1:10],'bo',label='RK2')
plt.xlabel('t')
plt.ylabel('y')
plt.legend(loc='upper right',frameon=False)
# fig1.savefig('/mnt/c/Users/rickl/Documents/NumericalAnalysis/RK2scalar.pdf',bbox_inches="tight")

plt.show()






