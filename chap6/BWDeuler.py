# -*- coding: utf-8 -*-
"""
Created on Tue May 22 18:30:35 2018

@author: luk
"""

import numpy as np
t0 = 0.0; t1 = 2.0  # 区间
y0 = 0.5  # 初值
h = 0.01  # 步长
n = round( (t1-t0)/h )
eps = 1e-6
nloop = 1000  # 最大循环次数

x = np.zeros( (n+1,), dtype = np.float64)
y = np.zeros( (n+1,), dtype = np.float64)
x[0] = t0
y[0] = y0

# func和DerivativeFunc依据初值问题而定
func = lambda h, y1, y2: 9.0 * h * y1**3 - 8.0 * h * y1**2 + ( 1.0 - h ) * y1 - y2
DerivativeFunc = lambda h, y1: 27.0 * h * y1**2 - 16.0 * h * y1 +  1.0 - h 

for i in range(1,n+1):
    j = 0
    y_tmp = y0  # 每一个节点初值的选取都是上一个节点的值，这样可以加快收敛速度
    while True:
        yy = y_tmp - func( h, y_tmp, y0 ) / DerivativeFunc( h, y_tmp )  # 使用牛顿迭代法求解每一个节点的值
        j = j + 1
        if ( abs(yy-y_tmp) < eps or j > nloop ): break
        y_tmp = yy
       
    y0 = yy
    y[i] = yy
    x[i] = t0 + i * h
    

# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, y, 'b--')

plt.title( "BackEuler Method" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["y"], loc = 2 )  # 设置图例, loc = 9 将legend的位置置于顶部中间

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()

print( "please input Enter and  stop!" )
input()  