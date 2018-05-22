# -*- coding: utf-8 -*-
"""
Created on Mon May 21 21:24:02 2018

@author: luk
"""

import numpy as np
m = 4; n = 3
# a = np.array( [ [1., 2., -2., 1.], [1., 3., -2., 1.], [2., 2., 6., 4.], [2., 2., -4., 5.], [-4., 5., 3., 3.] ], dtype = np.float64 )
a = np.array( [ [1., 3., -2.], [2., 2., 6.], [2., 2., -4.], [-4., 5., 3.] ], dtype = np.float64 )
# a = np.array( [ [1., 3.], [2., 2.], [2., 2.] ], dtype = np.float64 )
Q = np.zeros( (m,n), dtype = np.float64 )
R = np.zeros( (n,n), dtype = np.float64 )
A0 = np.zeros_like( a, dtype = np.float64 )
y = np.zeros( (m,), dtype = np.float64 )

for j in range(n):  
    y[:] = a[:,j]
    for i in range(j):  
        R[i,j] = np.dot(Q[:,i],y[:])
        y[:] = y[:] - R[i,j] * Q[:,i]
    
    R[j,j] = np.sqrt( np.dot(y,y) )
    Q[:,j] = y / R[j,j]
    
print( "The matrix Q is:" )
for i in range(m):
    print( "%9f  "*len(Q[i,:]) % tuple(Q[i,:]) )
    
print( "The matrix R is:" )
for i in range(n):
    print( "%9f  "*len(R[i,:]) % tuple(R[i,:]) )
    
A0 = np.matmul( Q, R )
print( "The matrix A0 is:" )
for i in range(m):
    print( "%9f  "*len(A0[i,:]) % tuple(A0[i,:]) )
    

# 调用python库求解QR
print( "------------------" )
print( "调用python库求解QR" )
QQ, RR = np.linalg.qr( a, mode = 'reduced' )
print( "The matrix Q is:" )
for i in range(m):
    print( "%9f  "*len(QQ[i,:]) % tuple(QQ[i,:]) )

print( "The matrix R is:" )
for i in range(n):
    print( "%9f  "*len(RR[i,:]) % tuple(RR[i,:]) )
    
AA = np.matmul( QQ, RR )
print( "The matrix AA is:" )
for i in range(m):
    print( "%9f  "*len(AA[i,:]) % tuple(AA[i,:]) )
    
print( "please input NUMBER and stop!" )
input()