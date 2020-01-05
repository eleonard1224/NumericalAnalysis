'''
Different quadrature routines based on Newton-Cotes (NC) interpolation
'''

# Importing the libraries
import math
import numpy as np

# Integration is to be performed on the function f over the interval [a,b]
# i.e. int_{a}^{b} f(x)dx
# If the integration is to be performed using m points including x1=a and xm=b, 
# then a polynomial interpolant of f is constructed of order m-1 using the m points
# and the integral I = int_{a}^{b} p_(m-1)(x)dx can be evaulated analytically

# NCWeights stores pre-defined weights w1, w2, ..., wm for m = 2,3,4,5
# so that I = (xm-x1)*[w1*f(x1) + w2*f(x2) + ... + wm*f(xm)] 
def NCWeights(m):
    
    if m==2:
        w = np.transpose(np.array([[1, 1]]))/2;
    if m==3:
        w = np.transpose(np.array([[1, 4, 1]]))/6;
    if m==4:
        w = np.transpose(np.array([[1, 3, 3, 1]]))/8;
    if m==5:
        w = np.transpose(np.array([[7, 32, 12, 32, 7]]))/90;
        
    return w;

# QNC computes the integral I by taking the inner product of NCWeights(m) with
# the evaluation of the function f at m evenly-spaced points x1=a, x2, ..., xm=b
def QNC(f,a,b,m):
    
    x = np.linspace(a,b,m);
    y = f(x);
    w = NCWeights(m);
    I = (b-a)*np.matmul(y,w);
    return I;

# CompQNC computes I by splitting the interval [a,b] into n segments
# and then splits each of those n segments into m sub-segments in which
# QNC(m) is performed for each sub-segment and the results summed for all
# the n segments to evaluate I
def CompQNC(f,a,b,m,n):
    
    x = np.linspace(a,b,m*n-(n-1));
    y = f(x);
    delta = (b-a)/n;
    w = NCWeights(m);
    I = 0;
    
    for i in range(1,n+1):
        i1 = (m-1)*(i-1); 
        i2 = i1 + m;
        I += np.matmul(y[i1:i2],w)
    I = I*delta;    
    return I;

# AdaptQNC calculates the integral of f on [a,b] using adaptive quadrature
# The estimate of the error variable is based on a comparison of the integral
# over the current region with 1 segment composed of m sub-segments versus a comparison
# of the integral over the current region with 2 segments composed of m sub-segments
def AdaptQNC(f,a,b,m,tol):
    
    d = m-((m+1)%2);
    A1 = CompQNC(f,a,b,m,1);
    A2 = CompQNC(f,a,b,m,2);
    error = (A1-A2)/(math.pow(2,d+1)-1);
    if abs(error) < tol:
        I = A2 + error;
    else:
        mid = 0.5*(a+b);
        I1 = AdaptQNC(f,a,mid,m,tol/2);
        I2 = AdaptQNC(f,b,mid,m,tol/2);
        I = I1+I2;
    return I;

# Testing the quadrature routines on I = int_{0}^{pi/2} sin x dx = 1
f = np.sin; a = 0; b = np.pi/2; m = 5; n = 4; tol = 0.01;
I = QNC(f,a,b,m);
print('QNC(I) = ' + str(I))
I = CompQNC(f,a,b,m,n);
print('CompQNC(I) = ' + str(I))
I = AdaptQNC(f,a,b,m,tol);
print('AdaptQNC(I) = ' + str(I))



    
 
