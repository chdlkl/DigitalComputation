# -*- coding: utf-8 -*-
"""
Created on Mon May 21 18:08:19 2018

@author: luk
"""

import numpy as np
n1 = 21  # 插值基节点
n2 = 51  # 待插值节点
t0= 0.0; t1 = 2.0  # 插值区间
x = np.zeros( (n1,), dtype = np.float64 )
y = np.zeros( (n1,), dtype = np.float64 )
x2 = np.zeros( (n2,), dtype = np.float64 )
y2 = np.zeros( (n2,), dtype = np.float64 )
y1 = np.zeros( (n2,), dtype = np.float64 )  # 调用内部函数计算

for i in range(n1):  # 构造插值基节点
    x[i] = np.float64(i) * ( t1-t0 ) / np.float64(n1-1)
    y[i] = np.float64(1.0) / ( np.float64(1.0) + np.float64(12.0) * x[i] * x[i] )
    
for i in range(n2):  # 待插值节点
    x2[i] = np.float64(i) * ( t1-t0 ) / np.float64(n2-1)

# --------------- -三次样条插值函数-----------------
def spln( n, t, m, n1, n2 ):  # 插值子程序
    global x, y, tmp
    f = np.zeros( (m,), dtype = np.float64 )
    f1 = np.zeros( (m,), dtype = np.float64 )  # 为f的导数
    
    s2 = np.zeros( (n,), dtype = np.float64 )
    e = np.zeros( (n2,), dtype = np.float64 )
    
    h = np.zeros( (n1,), dtype = np.float64 )
    dy = np.zeros( (n1,), dtype = np.float64 )
    s = np.zeros( (n,), dtype = np.float64 )
   
    for i in range(n1):
        h[i] = x[i+1] - x[i]
        dy[i] = ( y[i+1] - y[i] ) / h[i]
        
    s2[0] = np.float64(0.); s2[n-1] = np.float64(0.)
    
    for i in range(1,n1):
        s2[i] = np.float64(6.) * ( dy[i] - dy[i-1] )
   
    z = np.float64(0.5) / ( h[0] + h[1] )

    s[0] = -h[1] * z
    e[0] = s2[1] * z

    for i in range(1,n2):
        k = i - 1
        j = i + 1
        z = np.float64(1.) / ( np.float64(2.) * ( h[i]+h[j] ) + h[i]*s[k] )
        s[i] = -h[j] * z
        e[i] = ( s2[j]-h[i]*e[k] ) * z
    
    s2[n1-1] = e[n2-1]
    
    for i in range(n2,1,-1):
        k = i - 1
        s2[i] = s[k]*s2[i+1] + e[k]
        
    for i in range(n1):
        s[i] = ( s2[i+1] - s2[i] ) / h[i]
        
    i = 1
    k = 0
    for j in range(m):
        while True:
            if ( t > x[i] ):
                k = i
                i = i + 1
            else:
                break
        h1 = t - x[k]
        h2 = t - x[i]
        h3 = h1 * h2
        h4 = s2[k] + h1*s[k]
        z = ( s2[i] + s2[k] + h4 ) / np.float64(6.)
        f[j] = y[k] + h1*dy[k] + h3*z
        f1[j] = dy[k] + z*( h1+h2 ) + h3 * s[k] / np.float64(6.)
    tmp = f[0]      
    return tmp
# --------------- -三次样条插值函数-----------------
    
# ----------------调用子程序进行插值----------------  
for i in range(n2):
    tmp = np.float64(0.)
    spln( n1, x2[i], 1, n1-1, n1-2 )
    y2[i] = tmp
# ----------------调用子程序进行插值---------------- 

# 调用函数插值
from scipy import interpolate 
tck = interpolate.splrep(x, y)
y1 = interpolate.splev(x2, tck)
    

# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, y, 'r.')
axes.plot(x2, y1, 'g--')
axes.plot(x2, y2, 'b--')

plt.title( "function: 1 / ( 1 + 12 * x^2 )" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["y","y1","y2"], loc = 9 )  # 设置图例, loc = 9 将legend的位置置于顶部中间

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()

print( "please input Enter and  stop!" )
input()