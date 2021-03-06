# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 18:54:14 2018

@author: luk
"""
#//---------------原方程-------------
#//     x + y + z + w = 10
#//    2x + 3y + z + w = 15
#//    3x - y + 2z - w = 3
#//    4x + y -3z + 2w = 5
#//----------------------------------
import numpy as np
m = 4 
eps = 1e-12
a = np.array( [ [1.0, 1.0, 1.0, 1.0], [2.0, 3.0, 1.0, 1.0], [3.0, -1.0, 2.0, -1.0], [4.0, 1.0, -3.0, 2.0] ], dtype = float )
b = np.array( [ 10.0, 15.0, 3.0, 5.0 ], dtype = float )
x = np.zeros( (m,), dtype = float )

print( '消去前矩阵a为：' )
for i in range( len(a) ):
      print( a[i,:] )
print( 'please input enter:' )
input()

for j in range( len(a)-1 ):
       if ( np.abs(a[j,j]) < eps ):
             print( ' the pivot is zero! ' )
             quit()
       for i in range(j+1,len(a)):
             mult = a[i,j] / a[j,j]
             for k in range( j,len(a) ):
                   a[i,k] = a[i,k] - mult * a[j,k]
             b[i] = b[i] - mult * b[j]
             
print( '消去后矩阵a为：' )
for i in range( len(a) ):
      print( a[i,:] )
print( 'please input enter:' )
input()

for i in range( len(a)-1,-1,-1 ):
      for j in range( i, len(a) ):
            b[i] = b[i] - a[i,j] * x[j]
      x[i] = b[i] / a[i,i]
      
print( '原方程解为：' )
for i in range(len(x)):
      print( x[i] )
      
print( 'please input enter:' )
input()

print( 'python求解：' )
a = np.array( [ [1.0, 1.0, 1.0, 1.0], [2.0, 3.0, 1.0, 1.0], [3.0, -1.0, 2.0, -1.0], [4.0, 1.0, -3.0, 2.0] ], dtype = float )
b = np.array( [ 10.0, 15.0, 3.0, 5.0 ], dtype = float )
x = np.zeros( (4,), dtype = float )
a = np.linalg.inv(a)  # 求矩阵a的逆矩阵
x = np.matmul(a,b)
for i in range(len(x)):
      print( x[i] )
