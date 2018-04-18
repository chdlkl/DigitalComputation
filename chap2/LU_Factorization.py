# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:02:30 2018

@author: luk
"""

# LU分解
# Ax = b; A = LU 
# LUx = b; Ux = c; Lc = b
import numpy as np
m = 3 # 方程组的阶数
a = np.array( [ [1,2,-1],[2,1,-2],[-3,1,1] ], dtype = np.float64 )
b = np.array( [3,3,-6], dtype = np.float64 )
c = np.zeros( (m,), dtype = np.float64 );  x = np.zeros( (m,), dtype = np.float64 )
L = np.zeros( (m,m), dtype = np.float64 ); U = np.zeros( (m,m), dtype = np.float64 )

def GetLU( a ):
      eps = 1e-6  # 当主元小于eps时，视主元为0，程序停止计算
      mult = 0.
      global L, U
      print( '原始矩阵a：' )
      for i in range(m):
            print( a[i,:] )
            # 对下三角矩阵L的对角线元素赋值
            L[i,i] = 1.0
      # 计算上三角矩阵U和下三角矩阵L
      U[:,:] = a[:,:]  # 此处不能写做U = a, 这样两者的id相同，两者数据同步变化
      for j in range(m-1):
            if ( np.abs( U[j,j] ) < eps ):
                  print( 'The pivot is zero!' )
                  quit()
            for i in range(j+1,m):
                  mult = U[i,j] / U[j,j]
                  L[i,j] = mult
                  for k in range(j,m):
                        U[i,k] = U[i,k] - mult * U[j,k]
      return L, U

GetLU( a )  # 调用函数计算LU矩阵
# 得到下三角矩阵L
print( 'L is:' )
for i in range(m):
      print( L[i,:] )

# 得到上三角矩阵U
print( 'U is:' )
for i in range(m):
      print( U[i,:] )
      
# 检验LU分解是否正确：a = LU
tmp_a = np.zeros( (m,m), dtype = np.float64 )
tmp_a = np.matmul(L,U)
for i in range(m):
      print( tmp_a[i,:] )
      
def BackSubstitution( b, c ):  # 求解x
      global L, U, x
      
      # 求c: Lc = b
      for i in range(m):
            for j in range(i):
                  b[i] = b[i] - L[i,j] * c[j]
            c[i] = b[i] / L[i,i]
      
      # 求x: Ux = c
      for i in range(m-1,-1,-1):
            for j in range(i+1,m):
                  c[i] = c[i] - U[i,j]*x[j]
            x[i] = c[i] / U[i,i]
      return x

BackSubstitution( b, c )  # 调用函数计算x
print( 'x is:' )
for i in range(m):
      print( x[i] )

# 检验x计算是否正确
b = np.matmul( a,x )
print( 'b is:' )
for i in range(m):
      print( b[i] )
print( 'please input Enter:' )
input()