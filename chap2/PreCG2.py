# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:57:51 2018

@author: luk
"""

"""
在病态矩阵问题上，共轭梯度法比部分主元的高斯消去法还要差
这种情况下，可通过预条件得到缓解，主要是将问题转化为良态矩阵系统，再实施CG方法
此代码的预条件因子M= ( D + wL ) * Inv(D) * ( D + wU )，为高斯-赛德尔预条件因子
A = L + D + U，L为A的下三角，U为A的上三角，D为A的对角线元素
"""
import numpy as np
n = 3
a = np.array( [[4., -2., 2. ],[-2., 2., -4.],[2., -4., 11.]], dtype = np.float64 )
b = np.array( [6., -10., 27.], dtype = np.float64 )
alpha = 0.0; beta = 0.0
M = np.zeros_like( a, dtype = np.float64 )
x = np.zeros( (n,), dtype = np.float64 )
d = np.zeros_like( x, dtype = np.float64 )
r = np.zeros_like( x, dtype = np.float64 )
z = np.zeros_like( x, dtype = np.float64 )
rtmp = np.zeros_like( x, dtype = np.float64 )
ztmp = np.zeros_like( x, dtype = np.float64 )

def GetM( a, n ):
    global M
    L = np.zeros_like( a, dtype = np.float64 )
    D = np.zeros_like( L, dtype = np.float64 )
    U = np.zeros_like( L, dtype = np.float64 )
    w = np.float64(1.0)  #w为高斯-赛德尔预条件因子，w是0-2之间的一个实数
    
    for i in range(n):
        D[i,i] = a[i,i]
        D[i,i] = np.float64(1.0) / D[i,i]
        
    for i in range(n):
        for j in range(n):
            if ( i < j ):
                U[i,j] = a[i,j]
            if ( i > j ):
                L[i,j] = a[i,j]
    M = D + w * L
    M = np.matmul( M, D )
    M = np.matmul( M, D+w*U )
    
    return M

GetM( a, n )

r = b - np.matmul( a,x )  # r0 = b - Ax0
z = np.matmul( M,r )  # z0 = Inv(M)*r0
d = z  # d0 = z0

for i in range(n):
    if ( all(r<0.0) ):
        quit()
    
    rtmp = r
    ztmp = z
    tmp = np.matmul( np.transpose(rtmp),ztmp ) / np.matmul( np.matmul( np.transpose(d),a ), d )
    alpha = tmp
    x = x + alpha * d
    r = r - alpha * np.matmul( a,d )
    z = np.matmul( M,r )
    tmp = np.matmul( np.transpose(r),z ) / np.matmul( np.transpose(rtmp),ztmp )
    beta = tmp
    d = z + beta * d
    
print( " the x of PreCG1 is :" )
for i in range(n):
    print( x[i] )
    
print( "the x of python is :" )  
a = np.array( [[4., -2., 2. ],[-2., 2., -4.],[2., -4., 11.]], dtype = np.float64 )
b = np.array( [6., -10., 27.], dtype = np.float64 )
a = np.linalg.inv(a)  # 求矩阵a的逆矩阵
x = np.matmul(a,b)
for i in range(len(x)):
      print( x[i] )
      
print( " please input Enter and stop!" )
input()