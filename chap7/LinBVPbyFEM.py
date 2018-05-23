# -*- coding: utf-8 -*-
"""
Created on Wed May 23 10:07:23 2018

@author: luk
"""

"""
y'' = 4y
y(0) = 1
y(1) = e^2
解析解: y = exp(2t)
"""

import numpy as np
n = 99 # 待求点的个数, 加上已知两个边界点，总共101个点
ya = 1.0; yb = np.exp(2.0)  # 边值
t0 = 0.0; t1 = 1.0  # 计算区域
h = (t1-t0)  / (n+1)  # 步长

tmpy = np.zeros( (n,1), dtype = np.float64 )
y = np.zeros( (n+2,1), dtype = np.float64 )  # 数值解
x = np.zeros( (n+2,1), dtype = np.float64 )
z = np.zeros( (n+2,1), dtype = np.float64 )  # 解析解

# 计算稀疏矩阵a: 对于不同的边值问题，系数矩阵a的元素要重新计算
# 对于y'' = f(y)[一般来说,y''=f(y',y,t)]来说，代码中的a是通用的
alpha = 8.0*h/3.0 + 2.0/h
beta = 2.0*h/3.0 - 1.0/h
a = np.zeros( (n,n), dtype = np.float64 )
for i in range(n):
    a[i,i] = alpha
    
for i in range(n-1):
    a[i,i+1] = beta
    a[i+1,i] = beta
    
# 计算右端项
b = np.zeros( (n,1), dtype = np.float64 )
b[0,0] = -ya*beta
b[n-1,0] = -yb*beta

# 解方程
a = np.linalg.inv(a)
tmpy = np.matmul(a,b)
y[0] = ya; y[n+1] = yb
y[1:n+1] = tmpy[:]

# 计算解析解
for i in range(n+2):
    x[i] = t0 + np.float64(i) * h
    z[i] = np.exp(2.0*x[i])
    
# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, y, 'r')
axes.plot(x, z, 'b--')

plt.title( "LinBVPbyFEM" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["y","z"], loc = 9 )  # 设置图例, loc = 9 将legend的位置置于顶部中间

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()

# 误差图
fig, axes = plt.subplots()
axes.plot(x, abs(z-y)/z*100, 'b--')

plt.title( "Error" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("error")

plt.show()

print( "please input Enter and  stop!" )
input()
    
    

