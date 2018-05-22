# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:15:02 2018

@author: luk
"""

import numpy as np
m = 4; n = 3
A = np.array( [ [1., 3., 2.], [2., 2., 3.], [2., 2., 4.], [-4., 5., 5.] ] )
Q = np.zeros( (m,m), dtype = np.float64 )
R = np.zeros( (m,n), dtype = np.float64 )

def calHouseholderQR():
    global A, m, n, Q, R
    II = np.zeros( (m,m), dtype = np.float64 )
    
    R[:] = A[:]
    k = 0
    for i in range(n):
        length = m - k
        x = np.zeros( (length,), dtype = np.float64 )
        w = np.zeros_like( x, dtype = np.float64 )
        v = np.zeros( (length,1), dtype = np.float64 )
        
        P = np.zeros( (length,length), dtype = np.float64 )
        tmpH = np.zeros( (length,length), dtype = np.float64 )
        tmpI = np.zeros( (length,length), dtype = np.float64 )
        
        
        x[:] = R[i:m,i]
        
        norm_x = np.dot( x,x )
        w[0] = np.sqrt( norm_x )  # w = [ |x|2, 0, ..., 0 ]
        v[:,0] = w - x  # v = w - x
        
        P = np.matmul( v, np.transpose(v) ) / np.dot( v[:,0], v[:,0] )  # P = v*v'/(v'*v)
        
        for j in range(length):
            tmpI[j,j] = np.float64(1.0)
            
        tmpH = tmpI - np.float64(2.0) * P  # H = I - 2P
        
        if (i < 1):
            R = np.matmul( tmpH, A )
            Q[:] = tmpH[:]
        else:
            II[:,:] = 0.0
            for j in range(m):
                II[j,j] = np.float64(1.0)
                
            II[i:m,i:m] = tmpH[:,:]
            Q = np.matmul( Q, II )
            R = np.matmul( II, R )
            
        k = k + 1

# 调用子程序计算QR
calHouseholderQR()
print( "The matrix Q is:" )
for i in range(m):
  print( "%9f  "*len(Q[i,:]) % tuple(Q[i,:]) )

print( "The matrix R is:" )
for i in range(m):
  print( "%9f  "*len(R[i,:]) % tuple(R[i,:]) )

# A = QR
A0 = np.matmul( Q, R )
print( "The matrix A is:" )
for i in range(m):
  print( "%9f  "*len(A0[i,:]) % tuple(A0[i,:]) )
  
print( "please input NUMBER and continue!" )
input()

# ---------------------------------------------------
print( "python求解QR" )
(QQ, RR) = np.linalg.qr( A, mode = 'complete' )
for i in range(m):
  print( "%9f  "*len(QQ[i,:]) % tuple(QQ[i,:]) )

print( "The matrix R is:" )
for i in range(m):
  print( "%9f  "*len(RR[i,:]) % tuple(RR[i,:]) )
  
AA = np.matmul( QQ, RR )
print( "The matrix A is:" )
for i in range(m):
  print( "%9f  "*len(AA[i,:]) % tuple(AA[i,:]) )
  
print( "please input NUMBER and stop!" )
input()