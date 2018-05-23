# -*- coding: utf-8 -*-
"""
Created on Tue May 22 21:43:31 2018

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
n = 3  # 已知前四个解 n=3,(0,1,2,3)表示四阶四步法
num = round( (t1-t0)/h )
eps = 1e-9

a = np.zeros( (n+1,), dtype = np.float64 )
b = np.zeros( (n+1,), dtype = np.float64 )
x = np.zeros( (num+1,), dtype = np.float64 )
u1 = np.zeros( (num+1,), dtype = np.float64 )  # 外插值
u2 = np.zeros( (num+1,), dtype = np.float64 )  # 内插值
u3 = np.zeros( (num+1,), dtype = np.float64 )  # 解析解

print( "求解Adams内插公式系数值" )
a[0] = 1.0
for i in range(1,n+1):
    s = 0.0; t = 2.0
    for j in range(i,0,-1):
        s = s + a[j-1] / t
        t = t + 1.0
    
    a[i] = 0.0 - s  # 此处与外插法不一样，a(i)=0.-s;而外插法为a(i)=1.-s
print( a ) 
  
print( "验证Adams内插公式系数值" )
sum1 = 0.0
for i in range(n+1):
    sum1 = sum1 + a[n-i] / (1+i)
    
if ( abs(sum1-0.0) < eps ):  # sum1=0.0,说明Adams内插公式系数分解正确
    print( "Adams内插公式系数分解正确" )
else:
    print( "Adams内插公式系数分解不正确" )
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

tmp = np.loadtxt("waicha.dat")
x[:] = tmp[:,0]; u1[:] = tmp[:,1]
u2[0:n] = u1[0:n]  # 前3个数值由外插算法给出

# 由于本程序用的是3阶，所以内插法的前3个数值要用外插法给出（wai插法要给出前4个数值）
for i in range(n,num+1):
    t11 = u1[i]
    while True:
        # 公式中的-3*u2由于初值问题中的y'=-3*y
        t22 = u2[i-1] + h * ( b[0]*(-3.0*t11) + b[1]*(-3.0*u2[i-1]) + b[2]*(-3.0*u2[i-2]) + b[3]*(-3.0*u2[i-3]) )
        if ( abs(t22-t11) < eps ):
            break
        else:
            t11 = t22  # 不满足迭代要求，继续迭代
    u2[i] = t22
    
# 计算解析解
for i in range(num+1):
    u3[i] = np.exp( -3.0*x[i] )
    
# 绘图
from matplotlib import pyplot as plt
fig, axes = plt.subplots()
axes.plot(x, u2, 'r')
axes.plot(x, u3, 'b--')

plt.title( "AdamsMoulton Method" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["u2","u3"], loc = 9 )  # 设置图例, loc = 9 将legend的位置置于顶部中间

# 调整loc的大小可改变其位置，loc的范围[0-10]
plt.show()

# 绘制误差图
fig, axes = plt.subplots()
axes.plot(x, (u1-u3)/u3*100, 'r')   # 外插百分比误差
axes.plot(x, (u2-u3)/u3*100, 'b--')  # 内插误差

plt.title( "AdamsBashforth and AdamsMoulton Error" )  # 设置图标题
plt.xlabel("x")  # 设置坐标轴名称
plt.ylabel("y")

axes.legend( ["u1","u2"], loc = 2 )  # 设置图例, loc = 9 将legend的位置置于顶部中间
plt.show()

print( "please input Enter and  stop!" )
input()   