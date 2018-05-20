# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:57:51 2018

@author: luk
"""

"""
[ 3 -1 0  0 0 0.5 ] [x1]   [2.5]
[ -1 3 -1 0 0.5 0 ] [x2]   [1.5]
[ 0 -1 3 -1  0  0 ] [x3] = [1.0]
[ 0  0 -1  3 -1 0 ] [x4]   [1.0]
[ 0 0.5 0 -1 3 -1 ] [x5]   [1.5]
[ 0.5 0 0 0 -1 3  ] [x6]   [2.5]
"""

import numpy as np
n = 6
maxInteration = 50
eps = 1e-12
a = np.array( (n,n), dtype = np.float64 )
b = np.array( (n,), dtype = np.float64 )

c = np.loadtxt( "IterationData.txt" )
a = c[0:n,0:n]
b = c[n,:]

def Jacobi( a, b, n, maxInteration, eps ):  # Jacobi迭代求解
    
    InvD = np.zeros_like( a, dtype = np.float64 )
    L = np.zeros_like( a, dtype = np.float64 )
    U = np.zeros_like( a, dtype = np.float64 )
    
    x0 = np.zeros_like( b, dtype = np.float64 )
    x = np.zeros_like( b, dtype = np.float64 )
    tmp = np.zeros_like( b, dtype = np.float64 )
    
    """
    !// Jacobi
    X0 = 初始向量
    Xk+1 = InvD[ b - (L+U)Xk ], k = 0, 1, 2,...
    InvD为系数矩阵对角线元素的逆矩阵
    b为右端项
    L为系数矩阵的下三角部分，注意与LU分解中的L不同
    U为系数矩阵的上三角部分，注意与LU分解中的U不同
    Xk为前一次计算出来的结果
    """
   
    for i in range(n):
        InvD[i,i] = 1.0 / a[i,i]
        
    for i in range(n):
        for j in range(n):
            if ( i > j ):
                L[i,j] = a[i,j]
            if ( i < j ):
                U[i,j] = a[i,j]
    
    i = 1
    while True:
        tmp = b - np.matmul( L+U, x0 )
        x = np.matmul( InvD, tmp )
        i = i + 1
        err = max( abs(x-x0) )
        if ( i > maxInteration or err < eps ):
            break
        x0 = x
    
    print( "Jacobi solution:" )
    print( x )
    print( "The iterations of Jacobi is ", i ) 
    
def Gauss_Seidel( a, b, n, maxInteration, eps ):  # 高斯-赛德尔迭代求解
    InvD = np.zeros_like( a, dtype = np.float64 )
    L = np.zeros_like( a, dtype = np.float64 )
    U = np.zeros_like( a, dtype = np.float64 )
    
    x0 = np.zeros_like( b, dtype = np.float64 )
    x = np.zeros_like( b, dtype = np.float64 )
    tmp = np.zeros_like( b, dtype = np.float64 )
    
    for i in range(n):
        InvD[i,i] = 1.0 / a[i,i]
        
    for i in range(n):
        for j in range(n):
            if ( i > j ):
                L[i,j] = a[i,j]
            if ( i < j ):
                U[i,j] = a[i,j]
    
    i = 1
    while True:
        tmp = b - np.matmul( L+U, x0 )  # X0为式中的Xk
        x = np.matmul( InvD, tmp )  # 先计算出右侧的Xk+1，右侧的Xk+1由Jacobi计算得到
        tmp = b - np.matmul( U, x0 ) - np.matmul( L, x )
        x = np.matmul( InvD, tmp )  # 最后计算出左侧的Xk+1
        i = i + 1
        err = max( abs(x-x0) )
        if ( i > maxInteration or err < eps ):
            break
        x0 = x
    
    print( "Gauss_Seidel solution:" )
    print( x )
    print( "The iterations of Gauss_Seidel is ", i ) 
    
def SOR( a, b, n, maxInteration, eps ):  # 连续过松弛迭代求解
    InvLD = np.zeros_like( a, dtype = np.float64 )
    D = np.zeros_like( a, dtype = np.float64 )
    L = np.zeros_like( a, dtype = np.float64 )
    LD = np.zeros_like( a, dtype = np.float64 )  # LD = wL + D
    U = np.zeros_like( a, dtype = np.float64 )
    
    x0 = np.zeros_like( b, dtype = np.float64 )
    x = np.zeros_like( b, dtype = np.float64 )
    tmp = np.zeros_like( b, dtype = np.float64 )
    
    w = np.float64(1.1)  # w为松弛因子，w大于0时，加快收敛（过松弛），小于0时，减缓收敛
    """
    X0 = 初始向量
    Xk+1 = Inv( wL + D ) * [ (1-w)*D*Xk - w*U*Xk ] + w*Inv( w*L + D )*b, k = 0, 1, 2,...
    """
    for i in range(n):
        D[i,i] = a[i,i]
        
    for i in range(n):
        for j in range(n):
            if ( i > j ):
                L[i,j] = a[i,j]
            if ( i < j ):
                U[i,j] = a[i,j]
                
    # 计算wL + D的逆矩阵
    LD = w*L + D
    InvLD = np.linalg.inv(LD)
    
    i = 1
    while True:
        tmp = np.matmul( ( np.float64(1.0) - w ) * D, x0 ) - np.matmul( w*U, x0 )  # (1-w)*D*Xk - w*U*Xk
        x = np.matmul( InvLD, tmp ) + np.matmul( w*InvLD, b )  # Xk+1 = Inv( wL + D ) * [ (1-w)*D*Xk - w*U*Xk ] + w*Inv( w*L + D )*b
        i = i + 1
        err = max( abs(x-x0) )
        if ( i > maxInteration or err < eps ):
            break
        x0 = x
    
    print( "SOR solution:" )
    print( x )
    print( "The iterations of SOR is ", i ) 
    
Jacobi( a, b, n, maxInteration, eps )
Gauss_Seidel( a, b, n, maxInteration, eps )
SOR( a, b, n, maxInteration, eps )

print( "python solution:" )
a = np.linalg.inv(a)  # 求矩阵a的逆矩阵
x = np.matmul(a,b)
for i in range(len(x)):
      print( x[i] )
      
print( " please input Enter and stop!" )
input()