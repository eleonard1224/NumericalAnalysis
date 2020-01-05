'''
Vandermonde interpolation
'''

# Importing the libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Defining the functions

# Given n data points (x1,y1), (x2,y2), ... , (xn,yn)
# vandermonde returns a vector (a1,a2,...,an) such that
# p_(n-1)(x) = a0+ a1*x + a2*x^2 + ... + a_(n-1)*x^(n-1) 
def vandermonde(x,y,n):
    
    A = np.ones((n,n))
    Acol = np.ones((n,1))
    
    for i in range(1,n):
        Acol = np.multiply(Acol,x)
        A[:,i] = Acol[:,0]
    
    a = np.linalg.solve(A,y)
    return a;

# Given a vector of a's, interpolant evaluates
# p_(n-1)(x) = a0+ a1*x + a2*x^2 + ... + a_(n-1)*x^(n-1)
# for linspace(x1,x2,d) and returns an dx2 matrix with 
# x values in the first column and y values in the second column 
def interpolant(a,x1,x2,n,d):

    x = np.linspace(x1,x2,d)

    intrp = np.zeros((d,2))
    intrp[:,0] = x
    
    A = np.ones((d,n))
    
    for i in range(1,n):
        x = np.multiply(A[:,i-1],x)
        A[:,i] = x

    y = np.matmul(A,a)
    
    intrp[:,1] = y[:,0]
    return intrp;

# Setting the data
n = 3
x = np.array([[5],[2],[6]])
y = np.array([[2],[4],[3]])
a = vandermonde(x,y,n)

x1 = 0; x2 = 8; d = 100;
intrp = interpolant(a,x1,x2,n,d)
print(intrp)

# Plotting the data
mpl.rcParams.update({'font.size': 16})
mpl.rc('lines', linewidth=4)
# Plot the results of the interpolation along with the interpolation points
fig1 = plt.figure(1)
plt.plot(intrp[:,0],intrp[:,1],'b',label='Interpolant')
plt.plot(x,y,'ks',markersize=8.5,label='Data') 
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='upper right',frameon=False)
# fig1.savefig('/mnt/c/Users/rickl/Documents/NumericalAnalysis/vandermonde.pdf',bbox_inches="tight")
plt.show()
    

