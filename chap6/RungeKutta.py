# -*- coding: utf-8 -*-
"""
Created on Tue May 22 18:49:06 2018

@author: luk
"""

"""
!// y' = t*y + t^3
!// y(t) = 3*exp(t^2/2) - t^2 - 2
"""
import numpy as np
t0 = 0.0; t1 = 1.0  # 区间
u0 = 1.0  # 初值
h = 0.01  # 步长
n = round( (t1-t0)/h )
t = np.zeros( (n+1,), dtype = np.float64)
u = np.zeros( (n+1,), dtype = np.float64)
y = np.zeros( (n+1,), dtype = np.float64)

derfunc = lambda t, u: t * u + t**3  # 由初值问题而定

u[0] = u0; t[0] = t0; y[0] = u0
for i in range(1,n+1):
    k1 = derfunc( t[i-1], u[i-1] )
    k2 = derfunc( t[i-1]+h/2.0, u[i-1]+h*k1/2.0 )
    k3 = derfunc( t[i-1]+h/2.0, u[i-1]+h*k2/2.0 )
    k4 = derfunc( t[i-1]+h, u[i-1]+h*k3 )
    t[i] = t0 + np.float64(i)*h
    u[i] = u[i-1] + h / 6.0 * ( k1 + 2.0*k2 + 2.0*k3 + k4 )
    y[i] = 3.0 * np.exp( t[i]**2/2.0 ) - t[i]**2 - 2.0
    
# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(t, u, 'r')
axes.plot(t, y, 'b--')

plt.title( "RungeKutta Method" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["u","y"], loc = 9 )  # 设置图例, loc = 9 将legend的位置置于顶部中间

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()

print( "please input Enter and  stop!" )
input()  
    