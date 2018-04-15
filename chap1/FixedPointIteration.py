# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 20:10:20 2018

@author: luk
"""

# 使用不动点迭代方法求解方程 x^3 + x - 1 = 0 的零点
eps = 1e-8
x0 = 0.5 # 初始迭代值

def func( x ):
    value = ( 1.0 + 2.0 * x**3 ) / ( 1.0 + 3.0 * x**2 )
    return value

while True:
    xtemp = func( x0 )
    if ( abs( xtemp - x0 ) < eps ):
        break
    x0 = xtemp

print ( " The solution is ", xtemp )
input()
