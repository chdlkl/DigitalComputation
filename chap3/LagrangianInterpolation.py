# -*- coding: utf-8 -*-
"""
Created on Sun May 20 19:39:01 2018

@author: luk
"""

import numpy as np
n1 = 51 # 插值基节点数  这里要特别注意，当n2的值过大时，会出现“龙格现象”！
t0 = np.float64(-1.); t1 = np.float64(1.)  # 插值区间
x = np.zeros( (n1,), dtype = np.float64 )
y = np.zeros( (n1,), dtype = np.float64 )

n2 = 501  # 待插值节点数  
xx = np.zeros( (n2,), dtype = np.float64 )
yy = np.zeros( (n2,), dtype = np.float64 )

for i in range(n1):
    x[i] = t0 + np.float64(i) * ( t1 - t0 ) / np.float64( n1 - 1 )
    y[i] = np.cos( x[i]*np.pi )
    
    
def GetInterpolation( x, y, t0, t1 ):
    global n1, n2, xx, yy
    L1 = np.zeros( (n1,), dtype = np.float64 )
    L2 = np.zeros( (n1,), dtype = np.float64 )
    
    for k in range(n2):
        xx[k] = t0 + np.float64(k) * ( t1 - t0 ) / np.float64( n2 - 1 )
        for i in range(n1):
            L1[i] = np.float64(1.0)
            for j in range(n1):
                if ( j != i ):
                    L1[i] = ( x[i] - x[j] ) * L1[i]  # 求插值分母
            
            L2[i] = np.float64(1.0)
            for j in range(n1):
                if ( j != i ):
                    L2[i] = ( xx[k] - x[j] ) * L2[i]  # 求插值分子
                    
        tmp = np.float64(0.0)
        for i in range(n1):
            tmp = tmp + y[i] * L2[i] / L1[i]  # 求取插值项
            
        yy[k] = tmp
    return xx, yy

# 拉格朗日插值法
GetInterpolation( x, y, t0, t1 )

# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, y, 'b.')
axes.plot(xx, yy, 'r-')
plt.title( "function: cos" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["y = cos(x)"], loc = 2 )  # 设置图例
plt.show()

print( "please input Enter and  stop!" )
input()
