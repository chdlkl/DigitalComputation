# -*- coding: utf-8 -*-
"""
Created on Tue May 22 15:52:36 2018

@author: luk
"""

# purpose: Romberg integral for ln(x) between 1 and 2
import numpy as np
n = 4  # 步长减半的次数
a = 1.0; b = 2.0  # 积分上下限
r = np.zeros( (n,n), dtype = np.float64 )

func = lambda x: np.log( x )  # 被积函数

# ------------------------------此部分不需要改--------------------------------
r[0,0] = ( b - a ) * ( func(a) + func(b) ) / 2.0
for j in range(1,n):
    
    total = 0.0
    h = ( b - a ) / 2**j
    tmp = int( 2**(j-1) + 1 )
    
    for i in range(1,tmp):
        total = total + func( a+(2*i-1)*h )
    
    r[j,0] = r[j-1,0] / 2.0 + h * total
    
    for k in range(1,j+1):
        r[j,k] = ( 4.0**k * r[j,k-1] - r[j-1,k-1] ) / ( 4.0**k - 1.0 )
        
# ------------------------------此部分不需要改--------------------------------
# output r
for i in range(n):
    print( "%12g"*len(r[i,:]) % tuple(r[i,:]) )
    
print( "积分为：", r[n-1,n-1] )
