# -*- coding: utf-8 -*-
"""
Created on Tue May 22 15:15:05 2018

@author: luk
"""

import numpy as np
a = 0.0; b = 1.0  # 积分区间, 被积函数为f(x) = x*e^x
n = 100 # 区间等份
h = ( b - a ) / ( n * 1.0 )
s = 0.0  # 积分结果

x = np.zeros( (n+1,), dtype = np.float64 )

for i in range(n+1):
    x[i] = a + np.float64(i) * h
    
# 计算辛普森积分
for i in range(1,n+1):
    y0 = x[i-1] * np.exp( x[i-1] )
    tmp = ( x[i-1] + x[i] ) / 2.0
    y1 = tmp * np.exp( tmp )
    y2 = x[i] * np.exp( x[i] )
    s = s + h/2.0 * ( y0 + 4.0 * y1 + y2 ) / 3.0  # 辛普森积分公式
    
print( "辛普森积分结果:", s )

# 计算梯形积分
s = 0.0
for i in range(1,n+1):
    y0 = x[i-1] * np.exp( x[i-1] )
    y1 = x[i] * np.exp( x[i] )
    s = s + h * ( y0 + y1 ) / 2.0  # 梯形积分公式
    
print( "梯形积分结果:", s )
print( "精确积分结果:", 1.0 )

print( "please input NUMBER and stop!" )
input()