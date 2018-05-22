# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:51:18 2018

@author: luk
"""

import numpy as np
# Nparas为模型系数的个数，ndata为观测数据，Niter为最大迭代次数
# 模型为 y = ( 0.362*x ) / ( 0.556+x )
Nparas = 2; ndata = 7; Niters = 50
xobs = np.array( [ 0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740 ], dtype = np.float64 )
yobs = np.array( [ 0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317 ], dtype = np.float64 )

def GS():
    global xobs, yobs, Nparas, ndata, Niters
    
    aInit = 0.1  # 初始猜测值
    bInit = 0.1  # 初始猜测值
    eps = 1.e-8  # 精度控制
    
    # yest为估计值，d为估计值和实际值yobs之间的残差，dT为数组d的变形，delta为增量矩阵
    yest = np.zeros( (ndata,), dtype = np.float64 )
    d = np.zeros_like( yest, dtype = np.float64 )
    dT = np.zeros( (ndata,1), dtype = np.float64 )
    delta = np.zeros( (Nparas,1), dtype = np.float64 )
    
    # J为雅可比矩阵，JT为J的转置矩阵
    J = np.zeros( (ndata,Nparas), dtype = np.float64 )
    JT = np.zeros( (Nparas,ndata), dtype = np.float64 )
    
    # Inv_H为H的逆矩阵，H为海森矩阵
    H = np.zeros( (Nparas,Nparas), dtype = np.float64 )
    Inv_H = np.zeros( (Nparas,Nparas), dtype = np.float64 )
    
    # aest，best分别为反演参数
    aest = aInit; best = bInit
    
    for it in range(Niters):
        yest = aest * xobs / ( best + xobs )  # 根据当前aest，best及xobs，得到函数值yest
        d = yobs - yest  # 计算已知值yobs与yest的残差
        
        for i in range(ndata):  # 计算雅可比矩阵。dy/da = x / ( b + x )，dy/db = -a*x / ( b + x )**2
            J[i,0] = xobs[i] / ( best + xobs[i] )
            J[i,1] = -aest * xobs[i] / ( best + xobs[i] )**2
            
        JT = np.transpose(J)
        H = np.matmul( JT,J )  # 计算海森矩阵
        
        # 计算步长delta，并根据步长计算新的参数估计值
        Inv_H = np.linalg.inv( H )
        
        dT[:,0] = d[:]
        delta = np.matmul( np.matmul(Inv_H,JT), dT )  # 计算增量
        atmp = aest + delta[0,0]
        btmp = best + delta[1,0]
        
        if ( np.dot(delta[:,0],delta[:,0]) < eps ): 
            break
        
        aest = atmp
        best = btmp
    return aest, best

# 调用函数计算    
(a,b) = GS()

n = 101  # 利用反演出来的计算计算
x = np.zeros( (n,), dtype = np.float64 )
y = np.zeros( (n,), dtype = np.float64 )
for i in range(n):
    x[i] = 0.0 + np.float64(i) * np.float64(4.0) / np.float64(n-1)  # 计算区间[0,4]
    y[i] = a * x[i] / ( b + x[i] )
    

# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(xobs, yobs, 'r.')
axes.plot(x, y, 'b--')
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

plt.show()

print( "please input NUMBER and stop!" )
input()
