'''
Function for returning a floating point number given a base, mantissa, and exponent
'''

# Importing the libraries
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Defining the floatpt method which returns a floating point number given a
#     base, mantissa, and exponent
# beta - base
# m - mantissa - 0.b1 b2 ... bt
# t - mantissa length
# e - exponent
def floatpt(beta,m,t,e):
    
    r = m * math.pow(beta,e-t);
    return r;

# Tests
r1 = floatpt(10,123,3,2);
print(r1);

r2 = floatpt(2,123,3,2);
print(r2);