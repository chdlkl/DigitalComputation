# -*- coding: utf-8 -*-
"""
Created on Tue May 22 16:29:18 2018

@author: luk
"""

import numpy as np
n = 11 # 设置求取高斯点的个数
eps = 1e-12  # 精度控制
a = -1.0; b = 1.0  # 积分区间

def N_Legendre(x):  # 生成n阶勒让德多项式
    global n
    
    a = np.zeros( (n,), dtype = np.float64 )
    a[0] = x  # 1阶勒让德多项式
    a[1] = 1.5 * x * x - 0.5  # 2阶勒让德多项式
    for i in range(2,n):
        # 利用递推关系产生n阶勒让德多项式
        a[i] = ( np.float64(i+i+1)*x*a[i-1] - np.float64(i)*a[i-2] ) / np.float64(i+1) 
    
    return a[n-1]

# print( N_Legendre(1.1) )
def N1_Legendre(x): # 生成n-1阶勒让德多项式
    global n
    
    a = np.zeros( (n,), dtype = np.float64 )
    a[0] = x
    a[1] = 1.5 * x * x - 0.5
    for i in range(2,n):
        a[i] = ( np.float64(i+i+1)*x*a[i-1] - np.float64(i)*a[i-2] ) / np.float64(i+1) 
        
    return a[n-2]
   
# print( N1_Legendre(1.1) )
    
def DN_Legendre(x): # 生成n阶勒让德多项式的导数表达式
    global n
    
    a = np.zeros( (n,), dtype = np.float64 )
    a[0] = x
    a[1] = 1.5 * x * x - 0.5
    for i in range(2,n):
        a[i] = ( np.float64(i+i+1)*x*a[i-1] - np.float64(i)*a[i-2] ) / np.float64(i+1) 
    
    result = ( a[n-2] - x*a[n-1] ) * n / ( 1.0 - x*x )
    
    return result
# print( DN_Legendre(1.1) )
    
def NewtonIteration( a, b ):  # 牛顿法求解函数的解
    nloop = 2000
    # a,b是传递进来的划分好的有一个解存在的区间
    x = ( a + b ) / 2.0
    i = 0
    while True:
        xtmp = x - N_Legendre(x) / DN_Legendre(x)  #  X(i+1) = Xi - f(Xi) / f'(Xi)  i = 1,2,...N
        i = i + 1
        if ( abs(xtmp-x) < eps and i > nloop ): break
        x = xtmp
    return x

f_root = np.zeros( (n,), dtype = np.float64 )  # 高斯点
f_coeff = np.zeros( (n,), dtype = np.float64 )  # 高斯系数

def root_coeff():  # 计算N阶勒让德多项式的根与去做权重系数
    global f_root, f_coeff
    global a, b
    
    h = 1e-6  # 步长
    j = 0  # 控制循环变量的初值
    m = -1.e0 - h 
    nstep = round( (b-a) / h )  # 此处的2为积分上限减去下限，积分区间为[-1,1]
    for i in range(nstep):
        if ( N_Legendre(m) * N_Legendre(m+h) < 0 ):
            f_root[j] = NewtonIteration( m, m+h )
            # 利用公式计算高斯点的权重
            f_coeff[j] = 2.0 / ( np.float64(n) * N1_Legendre(f_root[j]) * DN_Legendre(f_root[j]) )
            print( "高斯点序号:", j )
            print( "高斯点:", f_root[j], ",", "高斯点权重:", f_coeff[j] )
            print( "------------------------------------------------------- " )
            j = j + 1
        
        m = m + h  # 计算完后向前计算

# 计算高斯点与高斯权系数           
root_coeff()

# -----------------被积函数------------------
func = lambda x: np.exp( -x*x/2.0 )

# ------------------求积分-------------------
ans = 0.0
# 一般区间[a,b]上的求积分公式
# Integral[f(x),a,b] = Integral[f( ((b-a)*t+b+a)/2 )] * ( b-a )/2. t为N阶勒让德多项式的根
for i in range(n):
    ans = ans + f_coeff[i] * func( (a+b)/2.0 + (b-a)*f_root[i]/2.0 )
    
ans = ans * ( b-a ) / 2.0
print( "高斯-勒让德求积分结果:", ans )

print("please input NUMBER and stop!")
input()
