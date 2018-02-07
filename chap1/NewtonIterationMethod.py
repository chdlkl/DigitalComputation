eps = 1e-8  # 误差限度
nloop = 1000  # 迭代上限

i = 0  # 迭代次数初始值
x0 = -0.7  # 初始估计值

def func( x0 ):
    value = x0**3 + x0 - 1.0
    return value

def DerivativeFunc ( x0 ):
    value = 3.0 * x0**2 + 1.0
    return value

while True:
    xtemp = x0 - func( x0 ) / DerivativeFunc( x0 )  # X(i+1) = X(i) - f(Xi)/f'(Xi) i = 1, 2, ...,N
    i = i + 1
    if ( abs(xtemp-x0) < eps or i > nloop ):
        break
    x0 = xtemp

print ( ' solution is ', xtemp )
    
