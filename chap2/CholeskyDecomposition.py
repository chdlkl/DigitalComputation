# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:57:51 2018

@author: luk
"""
"""
首先给出对称正定矩阵的概念
如果A'= A，则n*n矩阵A是对称矩阵。如果对于所有向量x/=0，x'Ax > 0，则称矩阵A是正定矩阵
对称正定矩阵是非奇异的
任意矩阵A的行列式为矩阵对应的特征值的乘积
对于对称正定矩阵，可使用楚列斯基分解方法，可分解为: A = R'R，其中R是一个上三角矩阵
使用楚列斯基分解，对于对称正定矩阵，和一般的矩阵相比，它们只有一半数量的独立元素，可以用一半的计算代价实现，并且仅仅使用一半的内存
"""

import numpy as np
n = 3
a = np.array( [[4., -2., 2. ],[-2., 2., -4.],[2., -4., 11.]], dtype = np.float64 )
b = np.array( [6., -10., 27.], dtype = np.float64 )
R = np.zeros( (n,n), dtype = np.float64 )
RR = np.zeros_like( R, dtype = np.float64 )
x = np.zeros_like( b, dtype = np.float64 )

def getR( a, n ):
    global R, RR
    
    for i in range(n):
        if ( a[i,i] < 0.0 ):
            print( "矩阵a不是正定对称矩阵，程序结束！" )
        R[i,i] = np.sqrt( a[i,i] )
        j = n - i - 1
        if ( j >= 0 ):
            u = np.zeros( (j,1), dtype = np.float64 )
            u[:,0] = a[i,i+1:] / R[i,i]
            R[i,i+1:] = u[:,0]
            a[i+1:,i+1:] = a[i+1:,i+1:] - np.matmul( u, np.transpose(u) )
    
    RR = np.transpose( R )
    return R, RR

getR( a, n )  # 将矩阵a分解为R和RR，a = RR*R

"""
Ax = b, A = R'R
R'Rx = b,令R'c = b, 求出c
最后用Rx = c求出x
"""
def getRoot( b, n ):
    global x
    global R, RR
    c = np.zeros( (n,), dtype = np.float64 )
    # R'c = b
    for i in range(n):
        for j in range(i):
            b[i] = b[i] - RR[i,j] * c[j]
        c[i] = b[i] / RR[i,i]
    print( c ) 
    # 求x: Rx = c
    for i in range( n-1,-1,-1 ):
        for j in range( i+1,n ):
            c[i] = c[i] - R[i,j] * x[j]
        x[i] = c[i] / R[i,i]
    
    return x

getRoot( b, n )

print( " the x of CholeskyDecomposition is: " )
for i in range(n):
    print( x[i] )

print( "the x of python is :" )  
a = np.array( [[4., -2., 2. ],[-2., 2., -4.],[2., -4., 11.]], dtype = np.float64 )
b = np.array( [6., -10., 27.], dtype = np.float64 )
a = np.linalg.inv(a)  # 求矩阵a的逆矩阵
x = np.matmul(a,b)
for i in range(len(x)):
      print( x[i] )
     
print( "please input Ente and stop!" )
input()