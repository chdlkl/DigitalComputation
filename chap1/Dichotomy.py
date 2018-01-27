# 二分法求解方程f(t) = cos(t) - t 在区间[0,1]上的解
n = 40;  # 迭代次数
a = 0.0; b = 1.0;  # 搜索区间

import math
fa = math.cos( a ) - a;
fb = math.cos( b ) - b;
if fa * fb > 0.0:
  print ( " error: f(a) * f(b) < 0 is not satisfied! " )

i = 1
while i < n:
  c = ( a + b ) / 2.0
  fc = math.cos( c ) - c
  if fc == 0.0:
    break
  if fa * fc > 0.0:
    a = c; fa = fc;
  elif fa * fc < 0.0:
    b = c; fb = fc;

t = ( b + a ) / 2.0
print ( " solution is ", t )
print ( " error is ", abs(math.cos(t) - t) )

