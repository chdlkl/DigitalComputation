x0 = 0.0; x1 = 1.0  # 迭代初始值
eps = 1e-8  # 误差限度
nloop = 1000  # 最大循环次数

def func( x ):
  result = x**3 + x - 1.0
  return result

i = 0  # 循环叠加器
while True:
  xtmp = x1 - func(x1) * ( x1 - x0 ) / ( func(x1) - func(x0) )
  # Xi+1 = Xi - f(Xi) * ( Xi - Xi-1 ) / ( f(Xi) - f(Xi-1) ) i = 1, 2, ...
  if ( i > nloop or abs( xtmp - x1 ) <= eps ):
    break
  x0 = x1
  x1 = xtmp
  i = i + 1

print ( " The root of x^3 + x - 1.0 is ", xtmp )
