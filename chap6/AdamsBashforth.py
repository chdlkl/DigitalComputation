# -*- coding: utf-8 -*-
"""
Created on Tue May 22 19:16:15 2018

@author: luk
"""

"""
y'=-3*y
y(0) = 1
t = [0,2]
解析解：y = e^(-3t)
"""

import numpy as np
t0 = 0.0; t1 = 2.0  # 求值区间
h = 0.01 # 步长
y0 = 1.0
n = 3  # 已知前四个解 n=3,(0,1,2,3)表示四阶四步法
num = round( (t1-t0)/h )
eps = 1e-6

a = np.zeros( (n+1,), dtype = np.float64 )
b = np.zeros( (n+1,), dtype = np.float64 )
x = np.zeros( (num+1,), dtype = np.float64 )
u1 = np.zeros( (num+1,), dtype = np.float64 )  
u2 = np.zeros( (num+1,), dtype = np.float64 )  

print( "求解Adams外插公式系数值" )
a[0] = 1.0
for i in range(1,n+1):
    s = 0.0; t = 2.0
    for j in range(i,0,-1):
        s = s + a[j-1] / t
        t = t + 1.0
    
    a[i] = 1.0 - s  
print( a ) 
  
print( "验证Adams外插公式系数值" )
sum1 = 0.0
for i in range(n+1):
    sum1 = sum1 + a[n-i] / (1+i)
    
if ( abs(sum1-1.0) < eps ):  # sum1=1.0,说明Adams外插公式系数分解正确
    print( "Adams外插公式系数分解正确" )
else:
    print( "Adams外插公式系数分解不正确" )
    quit()
  
def calculate_jL( j, L ):
    global jL
    
    if ( L == 0 ):
        jL = 1.0
        return
    else:
        product1 = 1.0; product2 = 1.0
        for i in range(1,L+1):
            product1 = product1 * np.float64(i)  # 求解分子s(s-1)...(s-j+1)
            
        for i in range(j,j-L,-1):
            product2 = product2 * np.float64(i)  # 求解j!
            
        jL = product2 / product1

    return
    
print( "求解系数Bkl值" )
for L in range(n+1):
    b1 = 0.0
    for j in range(L,n+1):
        jL = 0.0
        calculate_jL(j,L)  # 求解(s(s-1)...(s-j+1)/j!)
        b1 = b1 + (-1)**L * a[j] * jL
    b[L] = b1
    
print( b*24 )  # 乘以24只是为了显示为整数

# 系数求取完毕，接下来用梯形法求取前四个数
derfunc = lambda y: -3.0 * y  # 由初值问题而定

def Trapezoid( n, h, t0 ):
    global u1, u2, y0, x
    
    u1[0] = y0; u2[0] = y0
    for i in range(1,n+1):
        t = t0 + h* (i-1)
        y = y0 + h * ( derfunc( y0 ) + derfunc( y0 + h*derfunc( y0 ) ) ) / 2.0
        y0 = y
        u1[i] = y
        x[i] = t + h
        u2[i] = np.exp( -3.0*x[i] )

Trapezoid( n, h, t0 )  # 求取前四个值     
# AdamsForth公式求解
for i in range(n+1,num+1):
    x[i] = t0 + i * h
    u2[i] = np.exp( -3.0*x[i] )  # 解析解
    u1[i] = u1[i-1] + h * ( b[0]*(-3.0*u1[i-1]) + b[1]*(-3.0*u1[i-2]) + b[2]*(-3.0*u1[i-3]) + b[3]*(-3.0*u1[i-4]) )  # 数值解
    
    
# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, u1, 'r')
axes.plot(x, u2, 'b--')

plt.title( "AdamsBashforth Method" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["u1","u2"], loc = 9 )  # 设置图例, loc = 9 将legend的位置置于顶部中间

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()

# 绘制误差图
fig, axes = plt.subplots()
axes.plot(x, (u1-u2)/u2*100, 'r-.')   # 外插百分比误差

plt.title( "AdamsBashforth Error" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["AdamsBashforth Error"], loc = 2 )  # 设置图例, loc = 9 将legend的位置置于顶部中间
plt.show()

print( "please input Enter and  stop!" )
input()   
    
        
    