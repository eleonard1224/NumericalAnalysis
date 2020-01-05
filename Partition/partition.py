'''
Tests direct vs. recursive location of a point given a partition
'''

# Importing the libraries
import math
import numpy as np


# Defining direct location of a point xstar given a partion x
# where return_value=1 => x1 <= xstar <= x2, etc.
def direct(xstar,x,n): 
    
    counter = 0;
    if xstar == x[n-1]:
        counter += 1;
        print('direct counter = ' + str(counter));
        return n-1;
    else:    
        i = 1;
        while xstar >= x[i]:
            counter += 1;
            i += 1;
        print('direct counter = ' + str(counter));
        return i;


# Defining recursive location of a point xstar given a partition x
# where return_value=1 => x1 <= xstar <= x2, etc.
def recursive(xstar,x,n):
    
    counter = 0;
    if xstar == x[n-1]:
        counter += 1;
        print('recursive counter = ' + str(counter));
        return n-1;
    else:      
        left = 1; right = n-1;
        mid = math.floor(0.5*(left+right));
        while right > left+1: 
            counter += 1;
            if xstar > x[mid]:
                left = mid;
            else:
                right = mid;
            mid = math.floor(0.5*(left+right));
        print('recursive counter = ' + str(counter));
        return left+1;


# Defining a partition
x1 = 0.0; xn = 1.0; n = 516;
x = np.linspace(x1,xn,n);

# Testing the location functions
xstar = 0.925;
index_direct = direct(xstar,x,n);
print('index_direct = ' + str(index_direct))
index_recursive = recursive(xstar,x,n);
print('index_recursive = ' + str(index_recursive))











