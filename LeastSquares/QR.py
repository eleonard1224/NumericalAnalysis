'''
Solving the least squares problem through QR factorization
'''

# Importing the libraries
import math
import numpy as np

# rotate determines the parameters c and s such that when the rotation matrix
# G = [c s; -s c] multiplies a vector x = [x1; x2] the result is a vector whose
# second component is zero i.e. Gx = y where y = [y1; 0]
def rotate(x1,x2):
    
    if x2 == 0:
        c = 1;
        s = 0;
    else:
        if abs(x2) >= abs(x1):
            cotangent = x1/x2;
            s = 1/math.sqrt(1 + cotangent*cotangent);
            c = cotangent*s;
        else:
            tangent = x2/x1;
            c = 1/math.sqrt(1 + tangent*tangent);
            s = tangent*c;    
    
    return c, s;

# RotMat returns a matrix which defines a rotation in the i,j plane
def RotMat(n,i,j,c,s):
    
    G = np.identity(n);
    G[i,i] = c; G[j,j] = c;
    G[i,j] = s; G[j,i] = -s;
    return G;

# QRRot returns the breakdown A = QR where A is an mxn matrix (m>=n),
# Q is an orthogonal mxm matrix and R is an upper triangular mxn matrix
def QRRot(A):
    
    m = A.shape[0];
    n = A.shape[1];
    Q = np.identity(m);
    R = np.copy(A);
 
    for j in range(0,n):    
        for i in range(m-1,j,-1): 
            cs = rotate(R[i-1,j],R[i,j]); c = cs[0]; s = cs[1];
            G = RotMat(m,i-1,i,c,s);
            R = np.matmul(G,R);
            Q = np.matmul(Q,np.matrix.transpose(G));

#             # Bug in this commented-out implementation
#             R[[i-1,i],j:n] = np.matmul(np.array([[c,s],[-s,c]]),R[[i-1,i],j:n]); 
#             Q[:,[i-1,i]] = np.matmul(Q[:,[i-1,i]],np.array([[c,-s],[s,c]]));

    return Q, R;

# UTriSol solves the matrix equation Ax=b where A is an nxn upper triangular matrix
def UTriSol(A,b):
    
    n = A.shape[0];
    x = np.zeros((n,1));
    for i in range(n-1,-1,-1):
        x[i] = b[i];
        for j in range(i+1,n):
            x[i] = x[i] - A[i,j]*x[j];
        x[i] = x[i]/A[i,i];
        
    return x;


# LSq takes the mxn, rank n matrix A (m>=n) and finds x in Ax=b such that
# the 2-norm ||Ax-b||_{2} is minimized through the factorization of A=QR
# where Q is an orthogonal matrix and R is an upper triangular matrix 
# Lsq also evaluates ||Ax-b||_{2}
def LSq(A,b):
    
    m = A.shape[0];
    n = A.shape[1];
    Q,R = QRRot(A);
    QT = np.transpose(Q);
    QTA = np.matmul(QT,A);
    QTb = np.matmul(QT,b);
    x = UTriSol(QTA[0:n,0:n],QTb[0:n])
    if m==n:
        res=0;
    else:
        res = np.linalg.norm(b[n:m],2);
    
    return x, res;

    
# Example
A = np.array([[2,0],[-1,1],[0,2]]);
print('A =')
print(A)
b = np.array([[1],[0],[-1]]);
print('b =')
print(b)
x, res = LSq(A,b);
print('x =')
print(x)
print('res =')
print(res)
      
      


