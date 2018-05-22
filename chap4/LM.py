# -*- coding: utf-8 -*-
"""
Created on Tue May 22 13:30:25 2018

@author: luk
"""

import numpy as np
# Nparas为模型系数的个数，ndata为观测数据，Niter为最大迭代次数
# 模型为 y = 20.5*exp(-0.24*x)
Nparas = 2; ndata = 9; Niters = 50
xobs = np.array( [ 0.25,0.5,1.0,1.5,2.0,3.0,4.0,6.0,8.0 ], dtype = np.float64 )
yobs = np.array( [ 19.306,18.182,16.126,14.302,12.685,9.978,7.849,4.857,3.005 ], dtype = np.float64 )

def GS():
    global xobs, yobs, Nparas, ndata, Niters
    
    aInit = 0.0  # 初始猜测值
    bInit = 0.0  # 初始猜测值
    e = 0.0
    eps = 1.e-8  # 精度控制
    lamda = 0.01  # 阻尼系数初始值
    v = 10  # 调节因子
    
    # yest为估计值，d为估计值和实际值yobs之间的残差，dT为数组d的变形，delta为增量矩阵
    yest = np.zeros( (ndata,), dtype = np.float64 )
    d = np.zeros_like( yest, dtype = np.float64 )
    dT = np.zeros( (ndata,1), dtype = np.float64 )
    delta = np.zeros( (Nparas,1), dtype = np.float64 )
    eye = np.zeros( (Nparas,Nparas), dtype = np.float64 )  # 单位矩阵
    
    # J为雅可比矩阵，JT为J的转置矩阵
    J = np.zeros( (ndata,Nparas), dtype = np.float64 )
    JT = np.zeros( (Nparas,ndata), dtype = np.float64 )
    
    # H为海森矩阵
    H = np.zeros( (Nparas,Nparas), dtype = np.float64 )
    H_LM = np.zeros( (Nparas,Nparas), dtype = np.float64 )
    Inv_H_LM = np.zeros( (Nparas,Nparas), dtype = np.float64 )
    
    # 计算的新值与新误差
    yest_LM = np.zeros( (ndata,), dtype = np.float64 )
    d_LM = np.zeros( (ndata,), dtype = np.float64 )
    
    for i in range(Nparas):
        eye[i,i] = 1.0
        
    order = 1
    # aest，best分别为反演参数
    aest = aInit; best = bInit
        
    for it in range(Niters):
        if ( order == 1 ):
            yest = aest * np.exp( -best*xobs )  # 根据当前aest，best及xobs，得到函数值yest
            d = yobs - yest  # 计算已知值yobs与yest的残差
        
            for i in range(ndata):  # 计算雅可比矩阵。dy/da = x / ( b + x )，dy/db = -a*x / ( b + x )**2
                J[i,0] = np.exp( -best*xobs[i] )
                J[i,1] = -aest * xobs[i] * np.exp( -best*xobs[i] )
             
            JT = np.transpose(J)
            H = np.matmul( JT,J )  # 计算海森矩阵
        if ( it == 0 ): e = np.dot( d, d )
        H_LM = H + lamda * eye
        
        # 计算步长delta，并根据步长计算新的参数估计值
        Inv_H_LM = np.linalg.inv( H_LM )
        
        dT[:,0] = d[:]
        delta = np.matmul( np.matmul(Inv_H_LM,JT),dT )  # 计算增量
        atmp = aest + delta[0,0]
        btmp = best + delta[1,0]
        
        
        if ( np.dot(delta[:,0],delta[:,0]) < eps ): break
        yest_LM = atmp*np.exp( -btmp*xobs )
        d_LM = yobs - yest_LM
        e_LM = np.dot( d_LM,d_LM )

        # 根据误差，决定如何更新参数和阻尼系数
        # 迭代成功时将lamda减小，否则增大lamda
        if ( e_LM < e ):
            lamda = lamda / v
            aest = atmp
            best = btmp
            e = e_LM
            order = 1
        else:
            order = 0
            lamda = lamda * v

    return aest, best

# 调用函数计算    
(a,b) = GS()

n = 101  # 利用反演出来的计算计算
x = np.zeros( (n,), dtype = np.float64 )
y = np.zeros( (n,), dtype = np.float64 )
for i in range(n):
    x[i] = 0.0 + np.float64(i) * np.float64(10.0) / np.float64(n-1)  # 计算区间[0,10]
    y[i] = a * np.exp( -b*x[i] )
    

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