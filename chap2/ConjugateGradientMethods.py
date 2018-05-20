# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:57:51 2018

@author: luk
"""

import numpy as np
n = 3
a = np.array( [[4., -2., 2. ],[-2., 2., -4.],[2., -4., 11.]], dtype = np.float64 )
b = np.array( [6., -10., 27.], dtype = np.float64 )
alpha = 0.
beta = 0.
x = np.zeros( (n,), dtype = np.float64 )
d = np.zeros( (n,), dtype = np.float64 )
r = np.zeros_like( x, dtype = np.float64 )
rtmp = np.zeros_like( r, dtype = np.float64 )

d = b - np.matmul( a,x )
r = d

for i in range(n):
    if ( all( r<0.0 ) ):
        quit()
    rtmp = r
    tmp = np.matmul( np.transpose(r),r ) / np.matmul( np.matmul( np.transpose(d), a ), d )
    alpha = tmp
    x = x + alpha * d
    r = r - alpha * np.matmul( a,d )
    tmp = np.matmul( np.transpose(r),r ) / np.matmul( np.transpose(rtmp), rtmp ) 
    beta = tmp
    d = r + beta * d

print( "the x of ConjugateGradientMethods is :" )  
for i in range(n):
    print( x[i] )
    
print( "the x of python is :" )  
a = np.linalg.inv(a)  # 求矩阵a的逆矩阵
x = np.matmul(a,b)
for i in range(len(x)):
      print( x[i] )
     
print( "please input Ente and stop!" )
input()

    