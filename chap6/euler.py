# -*- coding: utf-8 -*-
"""
Created on Tue May 22 18:16:09 2018

@author: luk
"""

"""
初值问题如下
y' = t*y + t^3
y(0) = y0
t = [0,1]
解析解为：y(t) = 3*exp(t^2/2) - t^2 - 2
"""

import numpy as np
a = 0.0; b = 1.0  # 区间
n = 500  # 份数
y0 = 1.0  # 初值

h = (b-a) / n
x = np.zeros( (n+1,), dtype = np.float64)
y = np.zeros( (n+1,), dtype = np.float64)
xx = np.zeros( (n+1,), dtype = np.float64)
yy = np.zeros( (n+1,), dtype = np.float64)

derfunc = lambda t, y: t * y + t**3  # 依据初值问题而定
x[0] = a
y[0] = y0
# 欧拉法求解
for i in range(1,n+1):
    t = a + h * np.float64(i-1)
    y[i] = y0 + h * derfunc( t, y0 )
    x[i] = t + h
    y0 = y[i]

# 解析解
for i in range(n+1):
    xx[i] = a + h * np.float64(i)
    t = xx[i]
    yy[i] = 3.0 * np.exp( t**2/2.0 ) - t**2 - 2.0
    

# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, y, 'r')
axes.plot(xx, yy, 'b--')

plt.title( "Euler Method" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["y","yy"], loc = 9 )  # 设置图例, loc = 9 将legend的位置置于顶部中间

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()

print( "please input Enter and  stop!" )
input()
    


