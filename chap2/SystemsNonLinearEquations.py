# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:57:51 2018

@author: luk
"""
"""
方程组：-u^3 + v = 0
         u^2 + v^2 - 1 = 0
即：f1(u,v) = -u^3 + v  f2(u,v) = u^2 + v^2 - 1
多元牛顿法的求解适用于小型方程组，对于大型不适用，因为需要提前指导导数
多变量牛顿法
x0 = 初始向量
DF(Xk) * s = -F(Xk)
Xk+1 = Xk + s, k = 0, 1, 2, ...
对于多元牛顿法来说，当方程组的阶数较小时，求解DF与F相对容易
但当方程组的阶数过大时，求解DF与F相当繁琐。
"""

import numpy as np
m = 2
maxLoop = 50
eps = 1e-12

DF = np.zeros( (m,m), dtype = np.float64 )
F = np.zeros( (m,), dtype = np.float64 )
s = np.zeros( (m,), dtype = np.float64 )
x0 = np.zeros( (m,), dtype = np.float64 )
x = np.array( [1.0, 2.0], dtype =np.float64 )

def GetDF( m, x ):
    global DF
    coff = np.array( ( [[-3., 1.], [2., 2. ]] ), dtype = np.float64 )
    coffuv = np.zeros_like( DF, dtype = np.float64 )
    
    coffuv[0,0] = x[0]**2; coffuv[1,0] = x[0]
    coffuv[0,1] = 1.e0; coffuv[1,1] = x[1]
    DF = coff * coffuv
    
def GetF( m, x ):
    global F
    
    # f1 = -u^3 + v
    # f2 = u^2 + v^2 - 1
    F[0] = -( -x[0]**3 + x[1] )
    F[1] = -( x[0]**2 + x[1]**2 - 1.e0 )
    
x0[:] = x[:]
for i in range( maxLoop ):
    GetDF( m, x0 )
    GetF( m, x0 )
    DF = np.linalg.inv(DF)
    s = np.matmul(DF,F)
    x = x0 + s
    if ( max(abs(x-x0))  < eps):
        break
    x0 = x
    
for i in range(m):
    print( x[i] )

print( "please input Enter and stop!" )
input()