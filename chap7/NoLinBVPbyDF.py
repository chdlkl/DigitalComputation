# -*- coding: utf-8 -*-
"""
Created on Wed May 23 10:17:03 2018

@author: luk
"""

"""
y'' = y' + cosy
y(0) = 0
y(pi) = 1
"""

import numpy as np
n = 99 # 待求点的个数, 加上已知两个边界点，总共101个点
ya = 0.0; yb = 1.0  # 边值
t0 = 0.0; t1 = np.pi  # 计算区域
h = (t1-t0)  / (n+1)  # 步长
nloop = 20 # 迭代次数

tmpy = np.zeros( (n,1), dtype = np.float64 )
yy = np.zeros( (n,1), dtype = np.float64 )
y = np.zeros( (n+2,1), dtype = np.float64 )  # 数值解
x = np.zeros( (n+2,1), dtype = np.float64 )

def Inv( yy, n, h, ya, yb ):
    global tmpy
    
    # 计算系数矩阵a: 具体问题，系数矩阵与右端项具体求解
    a = np.zeros( (n,n), dtype = np.float64 )
    for i in range(n):
        a[i,i] = -2.0 + h*h*np.sin( yy[i]*np.pi/180.0 )
        
    for i in range(n-1):
        a[i,i+1] = 1.0 - h / 2.0
        a[i+1,i] = 1.0 + h / 2.0

    # 计算右端项
    b = np.zeros( (n,1), dtype = np.float64 )
    b[0,0] = ( 1.0 + h/2.0 ) * ya - 2.0 * yy[0] + ( 1 - h/2.0 ) * yy[1] - h*h*np.cos( yy[0]*np.pi/180.0 )
    b[n-1,0] = ( 1.0 + h/2.0 ) * yy[n-2] - 2.0 * yy[n-1] + ( 1 - h/2.0 ) * yb - h*h*np.cos( yy[n-1]*np.pi/180.0 )
    for i in range(1,n-1):
        b[i,0] = ( 1.0 + h/2.0 ) * yy[i-1] - 2.0 * yy[i] + ( 1 - h/2.0 ) * yy[i+1] - h*h*np.cos( yy[i]*np.pi/180.0 )
       
    
    # 解方程
    a = np.linalg.inv( a )
    tmpy = np.matmul(a,b)
    return tmpy
    
for i in range(nloop):
    yy[:] = y[1:n+1]
    Inv( yy, n, h, ya, yb )
    y[1:n+1] = y[1:n+1] - tmpy[:]  # 牛顿迭代法求解非线性方程组

y[0] = ya; y[n+1] = yb
for i in range(n+2):
    x[i] = t0 + np.float64(i) * h

# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, y, 'b-.')

plt.title( "NoLinBVPbyDF" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()  

print( "please input Enter and  stop!" )
input()
      
    