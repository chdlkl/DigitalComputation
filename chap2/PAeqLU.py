# -*- coding: utf-8 -*-
"""
Created on Tue May  1 21:30:38 2018

@author: luk
"""

# PA = LU，与部分主元的区别在于，求得一个置换矩阵P，作用于b
# 而在部分主元里，矩阵A和b是同时进行行变换的，PA = LU 中，A进行行变换，b的行变换由P来实现
#--------Ax = b-------
#-------PAx = Pb------
#--------PA = LU------
#-------LUx = Pb------
#---------求解--------
#-------Lc = Pb-------
#--------Ux = c-------
#---------解得x-------
import numpy as np
m = 3 # 方程组的阶数
a = np.array( [ [2,1,5],[4,4,-4],[1,3,1] ], dtype = np.float64 )
b = np.array( [5,0,6], dtype = np.float64 )
c = np.zeros( (m,), dtype = np.float64 );  x = np.zeros( (m,), dtype = np.float64 )
L = np.zeros( (m,m), dtype = np.float64 ); U = np.zeros( (m,m), dtype = np.float64 )
Pb = np.zeros( (m,), dtype = np.float64 ); 

def GetLU( a, b ):
    
      P = np.zeros( (m,m), dtype = np.float64 )
      tmpA = np.zeros( (m,m), dtype = np.float64 )
      arrP = np.zeros( (m,), dtype = np.float64 )
      arrA = np.zeros( (m,), dtype = np.float64 )
      
      tmpA[:,:] = a[:,:]
      global L, U, Pb
      # 构造置换矩阵P
      for i in range(len(P)):
          P[i,i] = 1.0
      
      
      # 计算置换矩阵P
      for j in range( len(P)-1 ): 
          arrP[:] = P[j,:]  # 切不可写成arrP = P. 矩阵赋值要特别留意
          arrA[:] = tmpA[j,:]
          tmp_P = tmpA[j:len(tmpA),j]
          re = np.where( tmp_P == np.max( np.abs(tmp_P) ) )
          kk = re[0]  # 将元组re转化为数组kk
          k = kk[0]  # 提取数组kk中的元素值
          k = k + j 
          P[j,:] = P[k,:]; tmpA[j,:] = tmpA[k,:]
          P[k,:] = arrP[:]; tmpA[k,:] = arrA[:]
      
      print( '置换矩阵P为：' )
      for i in range( len(P) ):
          print( P[i,:] )
      
      # 计算Pa与Pb
      Pa = np.matmul( P,a )
      Pb = np.matmul( P,b )
      
      print( '置换矩阵Pa：' )
      for i in range(m):
          print( Pa[i,:] )
          # 对下三角矩阵L的对角线元素赋值
          L[i,i] = 1.0
    
      # 计算上三角矩阵U和下三角矩阵L
      U[:,:] = Pa[:,:]  # 此处不能写做U = Pa, 这样两者的id相同，两者数据同步变化
      for j in range(m-1):
            for i in range(j+1,m):
                  mult = U[i,j] / U[j,j]
                  L[i,j] = mult
                  for k in range(j,m):
                        U[i,k] = U[i,k] - mult * U[j,k]
      return L, U, Pb

GetLU( a, b )  # 调用函数计算LU矩阵
# 得到下三角矩阵L
print( 'L is:' )
for i in range(m):
      print( L[i,:] )

# 得到上三角矩阵U
print( 'U is:' )
for i in range(m):
      print( U[i,:] )
      
# 检验LU分解是否正确：Pa = LU
print( '检验LU分解是否正确，Pa = LU' )
tmp_Pa = np.matmul(L,U)
for i in range(m):
      print( tmp_Pa[i,:] )
      
def BackSubstitution( Pb, c ):  # 求解x
      global L, U, x
      
      # 求c: Lc = Pb
      for i in range(m):
            for j in range(i):
                  Pb[i] = Pb[i] - L[i,j] * c[j]
            c[i] = Pb[i] / L[i,i]
      
      # 求x: Ux = c
      for i in range(m-1,-1,-1):
            for j in range(i+1,m):
                  c[i] = c[i] - U[i,j]*x[j]
            x[i] = c[i] / U[i,i]
      return x

BackSubstitution( Pb, c )  # 调用函数计算x
print( 'x is:' )
for i in range(m):
      print( x[i] )

# 检验x计算是否正确
b = np.matmul( a,x )
print( '利用求出来的x计算b，检验b是否正确，b is:' )
for i in range(m):
      print( b[i] )
print( 'please input Enter:' )
input()