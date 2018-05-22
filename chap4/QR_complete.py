# -*- coding: utf-8 -*-
"""
Created on Mon May 21 21:57:07 2018

@author: luk
"""

import numpy as np
m = 4; n = 2
# A = np.array( [[1.,3.,2.], [2.,2.,3.], [2.,2.,4.], [-4.,5.,5.]], dtype = np.float64 )
A = np.array( [[1.,3.], [2.,2.], [2.,2.], [-4.,5.]], dtype = np.float64 )
A0 = np.zeros( (m,n), dtype = np.float64 )
R = np.zeros( (m,n), dtype = np.float64 )
Q = np.zeros( (m,m), dtype = np.float64 )

R[:] = A[:]

def calQR():
  global Q, R, m, n
  I1 = np.zeros( (m,m), dtype = np.float64 )
  Q1 = np.zeros( (m,m), dtype = np.float64 )
  U = np.zeros( (m,1),  dtype = np.float64 )

  for i in range(m):
    Q[i,i] = np.float64(1.0)
    I1[i,i] = np.float64(1.0)

  for j in range(n):
    U[:,0] = np.float64(0.0)
    for i in range(j,m):
      U[i,0] = R[i,j]

    U[j,0] = U[j,0] + np.sqrt( np.dot(U[:,0],U[:,0]) ) * ( U[j,0] / abs(U[j,0]) )
    UTU = np.matmul( np.transpose(U), U )
    Q1 = I1 - 2.0 * np.matmul( U, np.transpose(U) ) / UTU
    R = np.matmul( Q1, R )
    Q = np.matmul( Q1, Q )
    Q = np.transpose( Q )

calQR()

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


