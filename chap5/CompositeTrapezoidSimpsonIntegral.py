# -*- coding: utf-8 -*-
"""
Created on Tue May 22 15:23:52 2018

@author: luk
"""

import numpy as np
a = 0.0; b = 1.0  # 积分区间, 被积函数为f(x) = x*e^x
n = 100 # 区间等份
h = ( b - a ) / ( n * 1.0 )
s = 0.0  # 积分结果

x = np.zeros( (n+1,), dtype = np.float64 )
y = np.zeros( (n-1,), dtype = np.float64 )
y2 = np.zeros( (n-1,), dtype = np.float64 )
y1 = np.zeros( (n,), dtype = np.float64 )

for i in range(n+1):
    x[i] = a + np.float64(i) * h
    if ( i > 0 and i < n ):
        y[i-1] = x[i] * np.exp( x[i] )
          
s = 0.0
# 计算复合梯形积分
y0 = x[0] * np.exp( x[0] )
ym = x[n] * np.exp( x[n] )
s = h  * ( y0 + ym + 2.0 * sum(y) ) / 2.0
print( "复合梯形积分积分结果:", s )

s = 0.0
# 计算复合辛普森积分
i = 0
for j in range(1,2*n+1):
    if ( np.mod( j,2 ) == 1 ):
        i = i + 1
        tmp = j * h / 2.0
        y1[i-1] = tmp * np.exp( tmp )
        
i = 0
for j in range(1,2*n):
    if ( np.mod(j,2) == 0 ):
        i = i + 1
        tmp = j * h / 2.0
        y2[i-1] = tmp * np.exp( tmp )
        
y0 = x[0] * np.exp( x[0] )
y2m = x[n] * np.exp( x[n] )

# 此时复合辛普森积分公式的步长是复合梯形积分的一半，所以为h/2.0
s = h / 2.0 * ( y0 + y2m + 4.0*sum(y1) + 2.0*sum(y2) ) / 3.0  # 复合辛普森积分公式
print( "复合辛普森积分积分结果:", s )    
print( "精确积分结果:", 1.0 )

print( "please input NUMBER and stop!" )
input()