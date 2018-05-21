# -*- coding: utf-8 -*-
"""
Created on Mon May 21 09:43:16 2018

@author: luk
"""

import numpy as np
n = 50  # 插值节点：牛顿插值法的插值基节点与待插值节点数一样多
a = np.float64(-3.0); b = np.float64(3.0)  # 插值区间
x = np.zeros( (n,), dtype = np.float64 )  # 插值基节点
y = np.zeros( (n), dtype = np.float64 )
c = np.zeros( (n), dtype = np.float64 )  # 存放系数
xx = np.zeros( (n,), dtype = np.float64 )  # 待插值节点
res = np.zeros( (n), dtype = np.float64 )  # 插值结果

for i in range(n):
    x[i] = a + np.float64(i) * ( b - a ) / np.float64( n - 1 )
    xx[i] = a / np.float(2.0) + np.float64(i) * np.float64( b - a ) / np.float(2.0) / np.float64( n - 1 )
    y[i] = np.exp( x[i] )

# -----------------子程序，不需修改------------------   
def newtdd( x, y ):  # 求解系数c
    global c, n
    v = np.zeros( (n,n), dtype = np.float64 )
    for i in range(n):
        v[i,0] = y[i]

    for j in range(1,n):
        for i in range(n-j):
            v[i,j] = ( v[i+1,j-1] - v[i,j-1] ) / ( x[i+j] - x[i] )

            
    c[:] = v[0,:]  # 得到系数
    return c
    
def CalInterpolation( x,xx,res ):  # 求解插值节点处的节点值
    global c, n
    d = np.zeros( (n,n), dtype = np.float64 )
    d[0,:] = np.float64(1.0)
    for j in range(n):
        for i in range(1,n):
            d[i,j] = d[i-1,j] * ( xx[j] - x[i-1] )
    
     
    for i in range(n):
        res[i] = np.matmul( c[:], d[:,i] )
    
    return res
    
newtdd( x, y )  # 求解系数c
CalInterpolation( x,xx,res )  # 求解插值节点处的节点值
# -----------------子程序，不需修改------------------   

# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, y, 'b.')
axes.plot(xx, res, 'r-')
plt.title( "function: Exp" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["y = exp(x)"], loc = 2 )  # 设置图例
plt.show()

print( "please input Enter and  stop!" )
input()